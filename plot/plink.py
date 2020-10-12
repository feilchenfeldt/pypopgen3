import matplotlib as mpl
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

class PlinkPCA(object):
    
    def __init__(self, filebase):
        self.filebase = filebase
        self._load_eigenvectors()
        self._load_eigenvalues()
    
    def _load_eigenvectors(self):
        eigenvectors = pd.read_csv(self.filebase + '.pca.eigenvec',
                                       sep=' ',index_col=0, header=None)
        eigenvectors.index.name = 'sample'
        eigenvectors.columns = ['sample'] + [i+1 for i in range(eigenvectors.shape[1]-1)]
        self.eigenvectors = eigenvectors
        
        
    def _load_eigenvalues(self):
        eigenvalues = []
        with open(self.filebase + '.pca.eigenval') as f:
            for line in f.readlines():
                eigenvalues.append(float(line.strip()))
        eigenvalues = np.array(eigenvalues)
        variance_explained = eigenvalues*100/np.sum(eigenvalues)
        self.eigenvalues = eigenvalues
        self.variance_explained = variance_explained

    @staticmethod    
    def add_legend(color_dic,**kwa):
        legend_elements = []
        for p,c in color_dic.items():
            legend_elements.append(mpl.lines.Line2D([0,0],[0,1],color=c,label=p, marker='o', linestyle='None',
                                          markersize=10, markeredgewidth=1.5))
            l = plt.legend(handles=legend_elements,**kwa)
        return l

    def plot_pca(self,pc_max=12, ncols=3,sample_to_color=None, color_dic=None):
        eigenvectors = self.eigenvectors
        
        nplots = pc_max // 2
        nrows = int(np.ceil(nplots/ncols))
        
        fig = plt.figure(figsize=(15,5*nrows))

        j=0
        
        if sample_to_color is not None:
            colors = [sample_to_color[s] for s in eigenvectors.index]
        else:
            colors = None
        
        
        
        for i in range(0, pc_max, 2):
            j+=1
            x = eigenvectors[i+1]
            y = eigenvectors[i+2]

            ax = fig.add_subplot(nrows,ncols,int((i+2)/2))
            ax.scatter(x,y,color=colors)
            legend_elements = []
            ax.set_xlabel("PC {} ({:.1f}% variance) ".format(i+1, self.variance_explained[i]))
            ax.set_ylabel("PC {} ({:.1f}% variance)".format(i+2, self.variance_explained[i+1]))
            if (i == pc_max-2) or (i == pc_max-1):
                if sample_to_color is not None:
                    legend = self.add_legend(color_dic)
                else:
                    legend = plt.legend()
        return ax, legend
