import pandas as pd
import numpy as np
import bioframe as bf

def MergeLoops(
    loopsDF, 
    cols1=["chrom1", "start1", "end1"], 
    cols2=["chrom2", "start2", "end2"], 
    iterations=3,
    fdrCol="FDR",
    sortLoops=True,
    randomSeed=1,
    flank=None
):
    """
        To merge overlapping loops. Loops with the lowest FDR will be selected
        
        NOTE: loops must be sorted and left anchor should come first than 
        right anchor
        
        PARAMETERS
        ----------
        loopsDF: pandas DataFrame
        cols1: list, default ["chrom1","start1","end1"]
            name of columns for left anchor
        cols2: list, default ["chrom2","start2","end2"]
            name of columns for left anchor
        iterations: int, default 3
            iteratively tries to merge the overlapping loops
        fdrCol: str, default "FDR"
            fdr column from mustache loop file
            if None then randomly select one
        sortLoops: bool, default True
            sort loops if not already sorted
        randomSeed: int, default 1
            random seed sampling
        flank: int, default None
            flank to add to both side of loops
        OUTPUT
        ------
        pandas DataFrame
        
    """
    df = loopsDF.copy()
    if fdrCol is None:
        fdrCol = "FDR"
        df[fdrCol] = np.random.rand(len(df))
    if flank is not None:
        df[cols1[1]] -= flank
        df[cols2[1]] -= flank
        df[cols1[2]] += flank
        df[cols2[2]] += flank
    if sortLoops:
        df = df.sort_values(by=[*cols1, 
                                fdrCol]
                           ).drop_duplicates(subset=[*cols1, 
                                                     *cols2]
                                            ).reset_index(drop=True)
    for itrCount in range(iterations):
        lAnchor = df[cols1].copy().rename(columns={cols1[0]:"chrom", 
                                                   cols1[1]:"start", 
                                                   cols1[2]:"end"}
                                         )
        rAnchor = df[cols2].copy().rename(columns={cols2[0]:"chrom", 
                                                   cols2[1]:"start", 
                                                   cols2[2]:"end"}
                                         )
        # randCol = np.random.randint(1,10,len(lAnchor))
        # lAnchor = lAnchor.assign(rand_col = randCol)
        # rAnchor = rAnchor.assign(rand_col = randCol)
        lAnchorOlp = bf.closest(lAnchor,
                                df2=None, 
                                return_index=True, 
                                # tie_breaking_col="rand_col", 
                                suffixes=('_ll', 
                                          '_lr'
                                         )
                               ).rename(columns={"distance":"distance_l"})
        rAnchorOlp = bf.closest(rAnchor, 
                                df2=None, 
                                return_index=True, 
                                # tie_breaking_col="rand_col", 
                                suffixes=('_rl', 
                                          '_rr'
                                         )
                               ).rename(columns={"distance":"distance_r"})
        merged = pd.concat((lAnchorOlp,rAnchorOlp),axis=1)
        overlappingLoops = merged.query("index_ll == index_rl &\
                                        index_lr == index_rr &\
                                        distance_l == 0 &\
                                        distance_r == 0", 
                                        engine="python"
                                       )
        olpIndexes = overlappingLoops[["index_ll","index_lr"]]
        # print(overlappingLoops.tail(10))#[["index_ll","index_lr","index_rl","index_rr"]]
        def _SwapOrder(a):
            if a[0] > a[1]:
                return([a[1],a[0]])
            else:
                return([a[0],a[1]])        
        indxsArray = np.unique(np.array(list(map(_SwapOrder, 
                                                 np.array(olpIndexes)
                                                )
                                            )
                                       ),axis=0)
        # olpIndexes = pd.DataFrame(indxsArray,columns=["index_ll","index_lr"])
        # print(indxsArray[:10])
        indxToDrop = []
        for i in indxsArray:
            # print(i[0],i[1],len(df))
            fdr1 = df.iloc[i[0]][fdrCol]
            fdr2 = df.iloc[i[1]][fdrCol]
            if fdr1 < fdr2:
                indxToDrop.append(i[1])
                # df.drop(index=i[1],inplace=True,errors="ignore")
            else:
                indxToDrop.append(i[0])
                # df.drop(index=i[0],inplace=True,errors="ignore")
        df.drop(index=np.unique(indxToDrop),inplace=True)
        df = df.reset_index(drop=True)
        df = df.sample(frac=1, random_state=randomSeed).reset_index(drop=True)
    df = df.sort_values(by=[*cols1,fdrCol]).reset_index(drop=True)
    
    if flank is not None:
        df[cols1[1]] += flank
        df[cols2[1]] += flank
        df[cols1[2]] -= flank
        df[cols2[2]] -= flank
    
    return df
    