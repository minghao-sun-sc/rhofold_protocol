
import torch
import torch.nn as nn
from torch import Tensor
from torch.nn import functional as F
from typing import Union, List

class MAE(nn.Module):
    def __init__(self, ):
        super().__init__()

    def forward(self, preds, targets):
        loss = F.l1_loss(preds, targets)
        return loss

class RMSD(nn.Module):
    def __init__(self) -> None:
        super(RMSD, self).__init__()

    def forward(self, coords_pred: List[Tensor], coords: List[Tensor]) -> Tensor:
        rmsds = []
        for coords_pred, coords in zip(coords_pred, coords):
            rmsds.append(torch.sqrt(torch.mean(torch.sum(((coords_pred - coords) ** 2), dim=1))))
        return torch.tensor(rmsds).mean()

class KabschRMSD(nn.Module):
    def __init__(self) -> None:
        super(KabschRMSD, self).__init__()

    def forward(self, coords_pred: List[Tensor], coords: List[Tensor], return_all = False, mode = 'global'):
        rmsds = []
        rmsds_local = []
        for coords_pred, coords in zip(coords_pred, coords):
            coords_pred_mean = coords_pred.mean(dim=0, keepdim=True)  # (1,3)
            coords_mean = coords.mean(dim=0, keepdim=True)  # (1,3)

            A = (coords_pred - coords_pred_mean).transpose(0, 1) @ (coords - coords_mean)

            U, S, Vt = torch.linalg.svd(A)

            corr_mat = torch.diag(torch.tensor([1, 1, torch.sign(torch.det(A))], device=coords_pred.device))
            rotation = (U @ corr_mat) @ Vt
            translation = coords_pred_mean - torch.t(rotation @ coords_mean.t())  # (1,3)

            coords = (rotation @ coords.t()).t() + translation

            rmsds.append(torch.sqrt(torch.mean(torch.sum(((coords_pred - coords) ** 2), dim=1))))
            rmsds_local.append(torch.sqrt(torch.sum(((coords_pred - coords) ** 2), dim=1)))

        if return_all:
            return rmsds if mode == 'global' else rmsds_local
        else:
            return torch.tensor(rmsds).mean() if mode == 'global' else rmsds_local


    # sel2: reference protein
    # """
    # # Select only C alphas
    # sel1 += ' and name CA'
    # sel2 += ' and name CA'
    # gdts = [] # Store GDT at the different thresholds
    # cutoffs =
    # cmd.align(sel1, sel2, cycles=0, transform=0, object='aln')
    # mapping = cmd.get_raw_alignment('aln') # e.g. [[('prot1', 2530), ('prot2', 2540)], ...]
    # distances = []
    # for mapping_ in mapping:
    #     atom1 = '%s and id %d'%(mapping_[0][0], mapping_[0][1])
    #     atom2 = '%s and id %d'%(mapping_[1][0], mapping_[1][1])
    #     dist = cmd.get_distance(atom1, atom2)
    #     cmd.alter(atom1, 'b = %.4f'%dist)
    #     distances.append(dist)
    # distances = numpy.asarray(distances)
    # gdts = []
    # for cutoff in cutoffs:
    #     gdts.append((distances <= cutoff).sum()/float(len(distances)))
    # out = numpy.asarray(zip(cutoffs, gdts)).flatten()
    # print ("GDT_%d: %.4f; GDT_%d: %.4f; GDT_%d: %.4f; GDT_%d: %.4f;"%tuple(out))
    # GDT = numpy.mean(gdts)
    # print ("GDT: %.4f"%GDT)
    # cmd.spectrum('b', 'green_yellow_red', sel1)
    # return GDT

class KabschGDT_TS(nn.Module):
    def __init__(self) -> None:
        super(KabschGDT_TS, self).__init__()

    def forward(self, coords_pred: List[Tensor], coords: List[Tensor]):
        rmsds = []

        for coords_pred, coords in zip(coords_pred, coords):
            coords_pred_mean = coords_pred.mean(dim=0, keepdim=True)  # (1,3)
            coords_mean = coords.mean(dim=0, keepdim=True)  # (1,3)

            A = (coords_pred - coords_pred_mean).transpose(0, 1) @ (coords - coords_mean)

            U, S, Vt = torch.linalg.svd(A)

            corr_mat = torch.diag(torch.tensor([1, 1, torch.sign(torch.det(A))], device=coords_pred.device))
            rotation = (U @ corr_mat) @ Vt
            translation = coords_pred_mean - torch.t(rotation @ coords_mean.t())  # (1,3)

            coords = (rotation @ coords.t()).t() + translation

            diff = torch.sqrt(torch.sum(((coords_pred - coords) ** 2), dim=1))

            r_list = []
            for t in [1., 2., 4., 8.]:
                r = torch.sum(diff < t).item() / diff.shape[0]
                r_list.append(r)
            rmsds.append(sum(r_list)/len(r_list))

        return rmsds

class RMSDmedian(nn.Module):
    def __init__(self) -> None:
        super(RMSDmedian, self).__init__()

    def forward(self, coords_pred: List[Tensor], coords: List[Tensor]) -> Tensor:
        rmsds = []
        for coords_pred, coords in zip(coords_pred, coords):
            rmsds.append(torch.sqrt(torch.mean(torch.sum(((coords_pred - coords) ** 2), dim=1))))
        return torch.median(torch.tensor(rmsds))


class RMSDfraction(nn.Module):
    def __init__(self, distance) -> None:
        super(RMSDfraction, self).__init__()
        self.distance = distance

    def forward(self, coords_pred: List[Tensor], coords: List[Tensor]) -> Tensor:
        rmsds = []
        for coords_pred, coords in zip(coords_pred, coords):
            rmsds.append(torch.sqrt(torch.mean(torch.sum(((coords_pred - coords) ** 2), dim=1))))
        count = torch.tensor(rmsds) < self.distance
        return 100 * count.sum() / len(count)


class CentroidDist(nn.Module):
    def __init__(self) -> None:
        super(CentroidDist, self).__init__()

    def forward(self, coords_pred: List[Tensor], coords: List[Tensor]) -> Tensor:
        distances = []
        for coords_pred, coords in zip(coords_pred, coords):
                distances.append(torch.linalg.norm(coords_pred.mean(dim=0)-coords.mean(dim=0)))
        return torch.tensor(distances).mean()


class CentroidDistMedian(nn.Module):
    def __init__(self) -> None:
        super(CentroidDistMedian, self).__init__()

    def forward(self, coords_pred: List[Tensor], coords: List[Tensor]) -> Tensor:
        distances = []
        for coords_pred, coords in zip(coords_pred, coords):
            distances.append(torch.linalg.norm(coords_pred.mean(dim=0)-coords.mean(dim=0)))
        return torch.median(torch.tensor(distances))


class CentroidDistFraction(nn.Module):
    def __init__(self, distance) -> None:
        super(CentroidDistFraction, self).__init__()
        self.distance = distance

    def forward(self, coords_pred: List[Tensor], coords: List[Tensor]) -> Tensor:
        distances = []
        for coords_pred, coords in zip(coords_pred, coords):
            distances.append(torch.linalg.norm(coords_pred.mean(dim=0)-coords.mean(dim=0)))
        count = torch.tensor(distances) < self.distance
        return 100 * count.sum() / len(count)


class MeanPredictorLoss(nn.Module):

    def __init__(self, loss_func) -> None:
        super(MeanPredictorLoss, self).__init__()
        self.loss_func = loss_func

    def forward(self, x1: Tensor, targets: Tensor) -> Tensor:
        return self.loss_func(torch.full_like(targets, targets.mean()), targets)


if __name__ == '__main__':
    # rmsd = RMSD()
    rmsd = KabschRMSD()
    cords = torch.randn([128,3])
    cords_gt = torch.randn([128,3])
    re = rmsd([cords], [cords_gt])
    print(re)
