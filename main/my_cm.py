import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def PrettyConfusionMatrix(cm, labels, save_path):
    
    cm_with_totals = np.zeros((cm.shape[0] + 1, cm.shape[1] + 1))
    cm_with_totals[:cm.shape[0], :cm.shape[1]] = cm
    cm_with_totals[-1, -1] = np.sum(cm)

    fig, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(cm_with_totals,  fmt='', cmap='Blues', cbar=False, cbar_kws={'drawedges': True})

    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j + 0.5, i + 0.5, f'{cm[i, j]}', color='black', ha='center', va='center', fontweight='bold')

    
    for i in range(cm.shape[0]):
        if i == 0:
            row_sum = np.sum(cm[i])
            precision = cm[i, i] / row_sum * 100
            ax.text(cm.shape[1] + 0.5, i + 0.5, f'{np.sum(cm[i, :])}', color='black', ha='center', va='center', fontweight='bold')    
            ax.text(cm.shape[1] + 0.5, i + 0.9, f'Precision: {precision:.1f}%', color='black', ha='center', va='center', fontsize=10)
        else:
            ax.text(cm.shape[1] + 0.5, i + 0.5, f'{np.sum(cm[i, :])}', color='black', ha='center', va='center', fontweight='bold')

    
    for j in range(cm.shape[1]):
        if j == 0:
            column_sum = np.sum(cm[:, j])
            recall = cm[j, j] / column_sum * 100
            ax.text(j + 0.5, cm.shape[0] + 0.5, f'{np.sum(cm[:, j])}', color='black', ha='center', va='center', fontweight='bold')  
            ax.text(j + 0.5, cm.shape[0] + 0.85, f'Recall/TPR/\nSensitivity: {recall:.1f}%', color='black', ha='center', va='center', fontsize=10)

        if j == 1:
            column_sum = np.sum(cm[:, j])
            specificity = cm[j, j] / column_sum * 100
            ax.text(j + 0.5, cm.shape[0] + 0.5, f'{np.sum(cm[:, j])}', color='black', ha='center', va='center', fontweight='bold')   
            ax.text(j + 0.5, cm.shape[0] + 0.85, f'TNR/\nSpecificity: {specificity:.1f}%', color='black', ha='center', va='center', fontsize=10)

    
    ax.text(cm.shape[1] + 0.5, cm.shape[0] + 0.5, f'{np.sum(cm)}', color='white', ha='center', va='center', fontweight='bold')

    f1_score = (2*precision*recall)/(precision+recall)
    ax.text(cm.shape[1] + 0.5, cm.shape[0] + 0.9, f'F1 Score: {f1_score:.1f}%', color='white', ha='center', va='center', fontsize=10)

    ax.xaxis.set_label_position('top')
    ax.set_xlabel('True Value', labelpad=15, ha='left', x=0.24, fontweight='bold', fontname='serif')
    ax.set_ylabel('Predicted Value', labelpad=15, y=0.65, fontweight='bold', fontname='serif')
    ax.xaxis.tick_top()
    ax.set_xticks(np.arange(cm.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(cm.shape[0]) + 0.5, minor=False)
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)

    plt.tight_layout()
    plt.savefig(save_path)
    plt.show()

