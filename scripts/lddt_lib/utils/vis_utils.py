from matplotlib import pyplot as plt
import numpy as np

def save_png(logit_dist,
             logit_c4=None,
             logit_o5=None,
             logit_n91=None,
             logit_p=None,
             save_path='tmp.png'):
    npz = logit_dist[:, :, 5:]
    DSTEP = 0.5
    dist_vals_np = (4.25 + DSTEP * np.arange(32, dtype=np.float32))
    dist_vals = dist_vals_np.reshape([1, 1, 32])
    preds_reg = np.sum(dist_vals * npz, axis=2)

    preds_reg[preds_reg <= 5] = 0
    preds_reg[preds_reg >= 19] = 0

    plt.figure(figsize=(12, 9))
    plt.subplot(3, 3, 1)
    plt.imshow(20 - preds_reg)
    plt.title('Predicted Distance')

    logit_tmp_ = np.argmax(logit_dist, 2)
    plt.subplot(3, 3, 2)
    plt.imshow(logit_tmp_)
    plt.title('C1 cls argmax')

    if logit_c4 is not None:
        logit_c4_tmp_ = np.argmax(logit_c4, 2)
        plt.subplot(3, 3, 4)
        plt.imshow(logit_c4_tmp_)
        plt.title('C4 cls argmax')

    if logit_o5 is not None:
        logit_o5_tmp_ = np.argmax(logit_o5, 2)
        plt.subplot(3, 3, 5)
        plt.imshow(logit_o5_tmp_)
        plt.title('O5 cls argmax')

    if logit_n91 is not None:
        logit_n91_tmp_ = np.argmax(logit_n91, 2)
        plt.subplot(3, 3, 6)
        plt.imshow(logit_n91_tmp_)
        plt.title('N91 cls argmax')

    if logit_p is not None:
        logit_n91_tmp_ = np.argmax(logit_p, 2)
        plt.subplot(3, 3, 7)
        plt.imshow(logit_n91_tmp_)
        plt.title('N91 cls argmax')

    plt.savefig(save_path)
    plt.close()
    print('    save png to: ' + save_path)