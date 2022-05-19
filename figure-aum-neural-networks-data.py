import numpy as np
import pandas as pd
from pytorch_lightning.metrics.classification import AUROC
import torch
import torchvision.datasets
compute_auc = AUROC()

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
            proportions=[("subtrain",0.8),("validation",0.2)]
            lengths=[int(p*len(data_set)) for s,p in proportions]
            sub_list=torch.utils.data.random_split(data_set, lengths)
            for (s,p),sub_set in zip(proportions,sub_list):
                data_dict[data_name][s] = sub_set
        else:
            data_dict[data_name]["test"]=data_set

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

def AUM(pred_tensor, label_tensor):
    """Area Under Min(FP,FN)

    Loss function for imbalanced binary classification
    problems. Minimizing AUM empirically results in maximizing Area
    Under the ROC Curve (AUC). Arguments: pred_tensor and label_tensor
    should both be 1d tensors (vectors of real-valued predictions and
    labels for each observation in the set/batch).

    """
    fn_diff = torch.where(label_tensor == 1, -1, 0)
    fp_diff = torch.where(label_tensor == 1, 0, 1)
    thresh_tensor = -pred_tensor.flatten()
    sorted_indices = torch.argsort(thresh_tensor)
    sorted_fp_cum = fp_diff[sorted_indices].cumsum(axis=0)
    sorted_fn_cum = -fn_diff[sorted_indices].flip(0).cumsum(axis=0).flip(0)
    sorted_thresh = thresh_tensor[sorted_indices]
    sorted_is_diff = sorted_thresh.diff() != 0
    sorted_fp_end = torch.cat([sorted_is_diff, torch.tensor([True])])
    sorted_fn_end = torch.cat([torch.tensor([True]), sorted_is_diff])
    uniq_thresh = sorted_thresh[sorted_fp_end]
    uniq_fp_after = sorted_fp_cum[sorted_fp_end]
    uniq_fn_before = sorted_fn_cum[sorted_fn_end]
    uniq_min = torch.minimum(uniq_fn_before[1:], uniq_fp_after[:-1])
    return torch.mean(uniq_min * uniq_thresh.diff())

loss_dict = {
    "logistic":torch.nn.BCEWithLogitsLoss(),
    "AUM":AUM,
    }

def compute_loss_pred(features, labels):
    pred_mat = model(features)
    pred_vec = pred_mat.reshape(len(pred_mat))
    return loss_fun(pred_vec, labels), pred_vec

loss_df_list = []
max_epochs=100
torch.manual_seed(1) #TODO multiple seeds.
# TODO for loop over loss functions.
loss_name = "AUM"
loss_fun = loss_dict[loss_name]
out_metrics = {
    "AUC":compute_auc,
    "loss":loss_fun
    }
# TODO for loop over learning rates.
# TODO for loop over data sets.
set_dict = data_dict["MNIST"]
model = LeNet5()
optimizer = torch.optim.Adam(model.parameters(), lr=0.0005)
for epoch in range(max_epochs):
    print(epoch)
    # first update weights.
    subtrain_loader = torch.utils.data.DataLoader(
        set_dict["subtrain"], shuffle=True, batch_size=50)
    for subtrain_features, subtrain_labels in subtrain_loader:
        optimizer.zero_grad()
        subtrain_loss, pred_vec = compute_loss_pred(
            subtrain_features, subtrain_labels)
        subtrain_loss.backward()
        optimizer.step()
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
                loss_df_list.append(pd.DataFrame({
                    "epoch":[epoch],
                    "loss_name":[loss_name],
                    "set_name":[set_name],
                    "out_name":[out_name],
                    "out_value":out_tensor.detach().numpy(),
                    }))
loss_df = pd.concat(loss_df_list)

