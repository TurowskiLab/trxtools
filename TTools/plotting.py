import matplotlib.pyplot as plt
import pandas as pd

def plotPCA(data=pd.DataFrame(), names=list(), title=str(), PClimit=1):
    nPCA = len([col for col in data.columns if 'PC' in col])
    axes = [ i +1 for i in range(nPCA)[:-1]]

    for nPC in axes:
        fig = plt.figure(figsize = (15 ,15))
        #         ax = fig.add_subplot(nPCA-1,1,nPC)
        ax = fig.add_subplot(1 ,1 ,1)
        a = data['PC ' +str(nPC)].tolist()
        b = data['PC ' +str(nPC+1)].tolist()

        ax.scatter(x=a ,y=b ,color='lightgray')

        if names:
            for i, txt in enumerate(data.index.tolist()):
                if txt in names:
                    ax.annotate(txt, (a[i], b[i]))

        #         for mark, c in zip(l1_fake,l1_fakeColor):
        #             markTemp = data.T[mark]
        #             ax.scatter(x=markTemp['PC'+str(nPC)],y=markTemp['PC'+str(nPC+1)],marker='X', color=c)

        ax.legend()
        ax.grid(True)
        #         ax.set_xlim(-1,1)
        plt.xlabel('PC ' +str(nPC))
        plt.ylabel('PC ' +str(nPC +1))
        plt.title(title)
        plt.show()

        if nPC==PClimit: break