# conda env torch2
import torch
torch.set_num_interop_threads(1)
torch.set_num_threads(1)
import pdb
import re
import pandas
import numpy
from sklearn.metrics import roc_auc_score
import sys

try:
    prog, seed_csv, loss_fun_name, step_size_str = sys.argv
except ValueError:
    print("Taking defaults!")
    seed_csv = "data_Classif_scaled/CIFAR10_N=10000_seed=1.csv"
    loss_fun_name = "SquaredHinge"
    step_size_str = '10000'
data_csv = re.sub("_seed=[0-9]", "", seed_csv)
step_size = float(step_size_str)
max_steps = 100000

Logistic = torch.nn.BCEWithLogitsLoss()
def SquaredHinge(predictions, labels, margin=1):
    torch.autograd.set_detect_anomaly(True)
    labels_length = len(labels)
    augmented_predictions = torch.where(labels == 1, predictions,
                                        predictions + margin)
    augmented_indices_sorted = torch.argsort(augmented_predictions)
    predicted_value = predictions[augmented_indices_sorted]
    labels_sorted = labels[augmented_indices_sorted]
    I_pos = torch.where(labels_sorted == 1)[0]
    I_neg = torch.where(labels_sorted == 0)[0]
    N = len(I_pos)*len(I_neg)
    z_coeff = margin - predicted_value
    a_coeff = torch.cumsum((labels_sorted)/N, dim = 0)
    b_coeff = torch.cumsum((2*z_coeff*labels_sorted)/N, dim = 0)
    c_coeff = torch.cumsum(((z_coeff**2)*labels_sorted)/N, dim = 0)
    loss_values = a_coeff*(predicted_value**2) + b_coeff*predicted_value + c_coeff
    return torch.mean(loss_values[I_neg])
def AUM(pred_tensor, label_tensor):
    """Area Under Min(FP,FN)

    Loss function for imbalanced binary classification
    problems. Minimizing AUM empirically results in maximizing Area
    Under the ROC Curve (AUC). Arguments: pred_tensor and label_tensor
    should both be 1d tensors (vectors of real-valued predictions and
    labels for each observation in the set/batch).

    """
    is_pos = label_tensor == 1
    is_neg = is_pos==False
    fp_diff = torch.where(is_neg, 1, 0)/is_neg.sum()
    fn_diff = torch.where(is_pos, -1, 0)/is_pos.sum()
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
    return torch.sum(uniq_min * uniq_thresh.diff())
loss_fun_names = ["Logistic","AUM","SquaredHinge"]
data_set_names = [
    #"STL10",
    "CIFAR10","MNIST","FashionMNIST"]
csv_names = [
    "data_Classif_scaled/CIFAR10_N=100.csv"
]
class LinearModel(torch.nn.Module):
    def __init__(self, ncol):
        super(LinearModel, self).__init__()
        self.weight_vec = torch.nn.Linear(ncol, 1, bias=False)
    def forward(self, feature_mat):
        return self.weight_vec(feature_mat)

weight_df=pandas.read_csv(seed_csv)
data_df = pandas.read_csv(data_csv)
y_tensor = torch.tensor(data_df.y).float()
X_tensor = torch.tensor(data_df.iloc[:,2:].to_numpy()).float()
X_nrow, X_ncol = X_tensor.shape
set_dict = {}
for set_name in "subtrain","validation":
    is_set = torch.tensor(data_df.set==set_name)
    set_dict[set_name] = {
        "X":X_tensor[is_set,],
        "y":y_tensor[is_set]
        }
loss_fun = eval(loss_fun_name)
model = LinearModel(X_ncol)
model.weight_vec.weight.data = torch.tensor(weight_df.weight).float()
optimizer = torch.optim.SGD(model.parameters(),lr=step_size)

out_df_list = []
for step_number in range(max_steps+1):
    loss_dict = {}
    for set_name, xy_dict in set_dict.items():
        set_nrow, set_ncol = xy_dict["X"].shape
        set_pred = model(xy_dict["X"]).reshape(set_nrow)
        set_loss = loss_fun(set_pred, xy_dict["y"])
        loss_dict[set_name] = set_loss
        set_auc = roc_auc_score(xy_dict["y"],set_pred.detach().numpy())
        out_row = pandas.DataFrame({
            'step_number':step_number,
            'set_name':set_name,
            "loss_value":set_loss.detach().numpy(),
            "auc":[set_auc],
        })
        print(out_row)
        out_df_list.append(out_row)
    if set_loss.isfinite()==False:
        break
    optimizer.zero_grad()
    loss_dict["subtrain"].backward()
    optimizer.step()
out_df = pandas.concat(out_df_list)
#"data_Classif_scaled/CIFAR10_N=100_seed=1.csv"
out_csv = seed_csv.replace("scaled","constant").replace(".csv", f"_lr={step_size}_loss={loss_fun_name}.csv")
out_df.to_csv(out_csv)
