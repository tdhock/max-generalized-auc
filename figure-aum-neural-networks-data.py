import pdb
import os
from pytorch_lightning.metrics.classification import AUROC
import torch
import torchvision.datasets
compute_auc = AUROC()
torch.set_num_interop_threads(1)#to avoid monsoon admin complaints.

class MySubset(torch.utils.data.Dataset):
    def __init__(self, dataset, indices):
        self.dataset = dataset
        self.indices = indices
    def __getitem__(self, index):
        i = self.indices[index].item()
        return self.dataset[i]
    def __len__(self):
        return len(self.indices)

data_dict = {}
data_name_tup = ("MNIST", "FashionMNIST")
set_tf_dict = {"train":True, "test":False}
for data_name in data_name_tup:
    data_class = getattr(torchvision.datasets, data_name)
    data_dict[data_name] = {}
    for in_name, train in set_tf_dict.items():
        data_set = data_class(
            ".", train=train, download=True,
            transform=torchvision.transforms.ToTensor(),
            target_transform=lambda label: 0 if label < 5 else 1)
        if in_name is "train":
            torch.manual_seed(1)
            label_list = [y for x,y in data_set]
            train_labels = torch.tensor(label_list)
            skip_dict = {0:1, 1:100}
            index_dict = {
                lab:(train_labels==lab).nonzero() for lab in skip_dict}
            index_list = [
                index_dict[lab][::skip] for lab,skip in skip_dict.items()]
            print([len(i) for i in index_list])
            indices = torch.cat(index_list)
            imbalanced = MySubset(data_set, indices)
            length_dict={"subtrain":int(0.8*len(imbalanced))}
            length_dict["validation"]=len(imbalanced)-length_dict["subtrain"]
            sub_list = torch.utils.data.random_split(
                imbalanced, length_dict.values())
            for s,sub_set in zip(length_dict.keys(),sub_list):
                data_dict[data_name][s] = sub_set
        else:
            data_dict[data_name]["test"]=data_set
dl = torch.utils.data.DataLoader(
    data_dict["MNIST"]["subtrain"], batch_size=2000)
for x,y in dl:
    pass
print(x.shape)
for data_name,set_dict in data_dict.items():
    for set_name,data_set in set_dict.items():
        print((data_name,set_name,len(data_set)))

class LeNet5(torch.nn.Module):
    def __init__(self):
        super(LeNet5, self).__init__()
        self.seq = torch.nn.Sequential(            
            torch.nn.Conv2d(
                in_channels=1, out_channels=6,
                kernel_size=5, stride=1, padding=2),
            torch.nn.ReLU(),
            torch.nn.AvgPool2d(kernel_size=2, stride=2),
            torch.nn.Conv2d(
                in_channels=6, out_channels=16,
                kernel_size=5, stride=1, padding=0),
            torch.nn.ReLU(),
            torch.nn.AvgPool2d(kernel_size=2, stride=2),
            torch.nn.Flatten(),
            torch.nn.Linear(in_features=400, out_features=120),
            torch.nn.ReLU(),
            torch.nn.Linear(in_features=120, out_features=84),
            torch.nn.ReLU(),
            torch.nn.Linear(in_features=84, out_features=1),
        )
    def forward(self, feature_mat):
        return self.seq(feature_mat)

def AUM(pred_tensor, label_tensor, rate=False):
    """Area Under Min(FP,FN)

    Loss function for imbalanced binary classification
    problems. Minimizing AUM empirically results in maximizing Area
    Under the ROC Curve (AUC). Arguments: pred_tensor and label_tensor
    should both be 1d tensors (vectors of real-valued predictions and
    labels for each observation in the set/batch).

    """
    is_positive = label_tensor == 1
    is_negative = label_tensor != 1
    fn_diff = torch.where(is_positive, -1, 0)
    fp_diff = torch.where(is_positive, 0, 1)
    thresh_tensor = -pred_tensor.flatten()
    sorted_indices = torch.argsort(thresh_tensor)
    fp_denom = torch.sum(is_negative) if rate else 1
    fn_denom = torch.sum(is_positive) if rate else 1
    sorted_fp_cum = fp_diff[
        sorted_indices].cumsum(axis=0)/fp_denom
    sorted_fn_cum = -fn_diff[
        sorted_indices].flip(0).cumsum(axis=0).flip(0)/fn_denom
    sorted_thresh = thresh_tensor[sorted_indices]
    sorted_is_diff = sorted_thresh.diff() != 0
    sorted_fp_end = torch.cat([sorted_is_diff, torch.tensor([True])])
    sorted_fn_end = torch.cat([torch.tensor([True]), sorted_is_diff])
    uniq_thresh = sorted_thresh[sorted_fp_end]
    uniq_fp_after = sorted_fp_cum[sorted_fp_end]
    uniq_fn_before = sorted_fn_cum[sorted_fn_end]
    uniq_min = torch.minimum(uniq_fn_before[1:], uniq_fp_after[:-1])
    return torch.sum(uniq_min * uniq_thresh.diff())

def AUM_rate(pred_tensor, label_tensor):
    """Area Under Min(FPR,FNR)"""
    return AUM(pred_tensor, label_tensor, rate=True)

out_cols = [
    "epoch",
    "step",
    "set_name",
    "out_name",
    "out_value"
    ]
loss_name = "logistic"
lr = 1e-3
seed=1
batch_size=50
def one_trial(loss_name, seed_str, lr_str, data_name, batch_size_str):
    out_csv = "/".join([
        "figure-aum-neural-networks-data",
        loss_name, seed_str, lr_str, data_name, batch_size_str,
        "steps.csv"
        ])
    out_dir = os.path.dirname(out_csv)
    os.makedirs(out_dir, exist_ok=True)
    def write_row(items,w_or_a):
        f=open(out_csv,w_or_a)
        f.write(",".join(items)+"\n")
    write_row(out_cols,"w")
    seed = int(seed_str)
    lr = float(lr_str)
    batch_size = int(batch_size_str)
    set_dict = data_dict[data_name]
    subtrain_label_list = [y for x,y in set_dict["subtrain"]]
    N_subtrain = len(subtrain_label_list)
    subtrain_label_tensor = torch.tensor(subtrain_label_list)
    label_count_dict = {
        lab:torch.sum(subtrain_label_tensor==lab) for lab in (0,1)}
    print(label_count_dict)
    label_weight_dict = {
        lab:N_subtrain/count for lab,count in label_count_dict.items()}
    def get_weight_tensor(lab_tensor):
        return torch.where(
            lab_tensor==0,
            label_weight_dict[0],
            label_weight_dict[1])
    torch.sum(get_weight_tensor(subtrain_label_tensor)) #should be 2.
    def log_loss(pred_tensor, label_tensor, *args):
        bce_inst = torch.nn.BCEWithLogitsLoss(*args)
        return bce_inst(pred_tensor, label_tensor.float())
    def weights_balanced(pred_tensor, label_tensor):
        return log_loss(
            pred_tensor, label_tensor,
            get_weight_tensor(label_tensor))
    loss_dict = {
        "logistic":log_loss,
        "balanced":weights_balanced,
        "AUM":AUM,
        "AUM_rate":AUM_rate,
    }
    loss_fun = loss_dict[loss_name]
    def compute_loss_pred(features, labels):
        pred_mat = model(features)
        pred_vec = pred_mat.reshape(len(pred_mat))
        return loss_fun(pred_vec, labels), pred_vec
    out_metrics = {
        "AUC":compute_auc,
        "loss":loss_fun
        }
    torch.manual_seed(seed) 
    model = LeNet5()
    optimizer = torch.optim.SGD(model.parameters(), lr=lr)
    for epoch in range(200):
        step = 0
        print(epoch)
        # first update weights.
        subtrain_loader = torch.utils.data.DataLoader(
            set_dict["subtrain"], shuffle=True, batch_size=batch_size)
        for subtrain_features, subtrain_labels in subtrain_loader:
            # then compute subtrain/validation loss.
            with torch.no_grad():
                for set_name,set_obj in set_dict.items():
                    set_loader = torch.utils.data.DataLoader(
                        set_obj, batch_size=100)
                    pred_list = []
                    label_list = []
                    for batch_features, batch_labels in set_loader:
                        batch_loss, batch_pred = compute_loss_pred(
                            batch_features, batch_labels)
                        pred_list.append(batch_pred)
                        label_list.append(batch_labels)
                    set_pred_tensor = torch.cat(pred_list)
                    set_label_tensor = torch.cat(label_list)
                    for out_name,out_fun in out_metrics.items():
                        out_tensor = out_fun(set_pred_tensor, set_label_tensor)
                        out_dict = {
                            "epoch":epoch,
                            "step":step,
                            "set_name":set_name,
                            "out_name":out_name,
                            "out_value":out_tensor.item()
                            }
                        item_list = [str(out_dict[N]) for N in out_cols]
                        write_row(item_list,"a")
            step += 1
            optimizer.zero_grad()
            subtrain_loss, pred_vec = compute_loss_pred(
                subtrain_features, subtrain_labels)
            subtrain_loss.backward()
            optimizer.step()

if __name__ == '__main__':
    import sys
    args = sys.argv[1:]
    print(args)
    one_trial(*args)
