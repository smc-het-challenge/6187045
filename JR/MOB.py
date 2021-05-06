## Modules to load
import pandas as pd
import numpy as np
from sklearn import mixture
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist,pdist
from sys import argv

# import data
df = pd.read_csv(argv[1],sep="\t")

# Consider only if #reads supporting mutation >10 and base calling with PHRED score >=30
#df = df[(df.tum_alt_allele_depth>10) & (df.tum_base_quality>=30)]

# add a colum for caryotype
df['caryotype'] = zip(df.nMin1_A,df.nMaj1_A)

####
##
## Choose the best gmm to model data
##
####
n_components_range = range(1,10)
cv_types = ['spherical', 'tied', 'diag', 'full']

def uniqify_list_ordered(seq): 
   # order preserving
   list_ordered = []
   for e in seq:
       if e not in list_ordered:
           list_ordered.append(e)
   return list_ordered


def selectGMM_chro(df,n_components_range,cv_types):
	gmm_clust = []
	chro_list = uniqify_list_ordered(df.chr)
	for chro in chro_list:
		sub_df = df[df.chr==chro]
		# select best model
		X = np.column_stack((sub_df.tum_allele_fraction, [1]* len(sub_df)))
		if n_components_range > len(X):
			n_components_range = range(1,len(X))
		lowest_bic = np.infty
		bic = []
		for cv_type in cv_types:
		    for n_components in n_components_range:
		        # Fit a mixture of Gaussians with EM
		        gmm = mixture.GMM(n_components=n_components,covariance_type=cv_type, verbose=0)
		        # Estimate model parameters for this EM algorithm
		        gmm.fit(X)
		        bic.append(gmm.bic(X))
		        # Choose this model if it abs(BIC) is the lowest 
		        if bic[-1] < lowest_bic:
		            lowest_bic = bic[-1]
		            best_gmm = gmm
		# predict data
		gmm = best_gmm
		gmm_clust = gmm_clust + list(gmm.fit_predict(sub_df[["tum_allele_fraction"]]))
	return gmm_clust



# with the selectGMM_chro function
df['gmm_clust'] = selectGMM_chro(df,n_components_range,cv_types)

# compute different possibillity of cluster proportion
df['clonalC1'] = 2*df.tum_allele_fraction / (1+2*df.tum_allele_fraction-(df.nMin1_A+df.nMaj1_A)*df.tum_allele_fraction)
df['clonalCMax'] = 2*df.tum_allele_fraction / (df.nMaj1_A+2*df.tum_allele_fraction-(df.nMin1_A+df.nMaj1_A)*df.tum_allele_fraction)
dfGpByClust = pd.DataFrame({
	'meanAF' : df.groupby(["chr","gmm_clust","nMin1_A","nMaj1_A","caryotype"])["tum_allele_fraction"].mean(),
	'Pclonal_Min' : df.groupby(["chr","gmm_clust","nMin1_A","nMaj1_A","caryotype"])["clonalC1"].mean(),
	'Pclonal_Maj' : df.groupby(["chr","gmm_clust","nMin1_A","nMaj1_A","caryotype"])["clonalCMax"].mean(),
	'mut_count' : df.groupby(["chr","gmm_clust","nMin1_A","nMaj1_A","caryotype"])["clonalCMax"].size()
	}).reset_index()

# do a kmean on the estimated proportion grouped by cluster
# first, construct a vector of estimated proportion (weighted by # of mutation)
forKmean = []
for i in range(len(dfGpByClust)):
	line = dfGpByClust.iloc[i]
	if line.Pclonal_Maj > 0 and line.Pclonal_Maj < 1:
		forKmean += [line.Pclonal_Maj]* line.mut_count
	if line.Pclonal_Min > 0 and line.Pclonal_Min < 1:
		forKmean += [line.Pclonal_Min]* line.mut_count

forKmean = np.column_stack((forKmean, [1]* len(forKmean)))

# choose the best k
# k range
k_range = range(2,max(dfGpByClust.gmm_clust)+1+1)
# fit kmeans model for each n_clusters = k
k_means_var = [KMeans(n_clusters=k).fit(forKmean) for k in k_range]
# pull out the cluster centers for each model 
centroids = [X.cluster_centers_ for X in k_means_var]
# calculate the Euclidean distance for each pointto each cluster center
k_euclid = [cdist(forKmean, cent, 'euclidean') for cent in centroids]
dist = [np.min(ke,axis=1) for ke in k_euclid]
#Total within-cluster sum of squares
wcss = [sum(d**2) for d in dist]
#total sum of square
tss = sum(pdist(forKmean)**2)/forKmean.shape[0]
# the between-cluster sum of squares
bss = tss - wcss

best_bss = 'inf'
for k in k_range:
	if bss[k-2] < best_bss:
		best_bss = bss[k-2]
		best_k = k

# do kmean with with the best k
kmean = KMeans(n_clusters=best_k)
kmean.fit(forKmean)

Pclonal_chosen = [1] *len(dfGpByClust)
K_chosen = [1] *len(dfGpByClust)
n_chosen = [1] *len(dfGpByClust)

# choose the n for each gmm cluster
for i in range(len(dfGpByClust)):
	line = dfGpByClust.iloc[i]
	if line.nMin1_A == line.nMaj1_A and line.nMin1_A == 1:
		n_chosen[i] = 1
		Pclonal_chosen[i] = line.Pclonal_Maj 
		K_chosen[i] = kmean.predict(np.array([line.Pclonal_Maj,1]).reshape(1, -1))[0]
	else :
		p_maj_dist = p_min_dist = ""
		if line.Pclonal_Maj > 0 and line.Pclonal_Maj < 1:
			p_maj_dist = min(min(kmean.transform(np.array([line.Pclonal_Maj,1]).reshape(1, -1))))
		if line.Pclonal_Min > 0 and line.Pclonal_Min < 1:
			p_min_dist = min(min(kmean.transform(np.array([line.Pclonal_Min,1]).reshape(1, -1))))
		if p_maj_dist and p_min_dist :
			if p_maj_dist <= p_min_dist :
				n_chosen[i] = dfGpByClust.iloc[i].nMaj1_A
				Pclonal_chosen[i] = dfGpByClust.iloc[i].Pclonal_Maj
				K_chosen[i] = kmean.predict(np.array([line.Pclonal_Maj,1]).reshape(1, -1))[0]
			else:
				n_chosen[i] = dfGpByClust.iloc[i].nMin1_A
				Pclonal_chosen[i] = dfGpByClust.iloc[i].Pclonal_Min
				K_chosen[i] = kmean.predict(np.array([line.Pclonal_Min,1]).reshape(1, -1))[0]
		elif p_maj_dist :
			n_chosen[i] = dfGpByClust.iloc[i].nMaj1_A
			Pclonal_chosen[i] = dfGpByClust.iloc[i].Pclonal_Maj
			K_chosen[i] = kmean.predict(np.array([line.Pclonal_Maj,1]).reshape(1, -1))[0]
		else :
			n_chosen[i] = dfGpByClust.iloc[i].nMin1_A
			Pclonal_chosen[i] = dfGpByClust.iloc[i].Pclonal_Min
			K_chosen[i] = kmean.predict(np.array([line.Pclonal_Min,1]).reshape(1, -1))[0]

dfGpByClust["Pclonal_chosen"] = Pclonal_chosen
dfGpByClust["K_chosen"] = K_chosen
dfGpByClust["K_chosen"] += 1
dfGpByClust["n_chosen"] = n_chosen

chro_list = uniqify_list_ordered(df.chr)

# to assign each cluster to each mutation
titi = ['']*len(df)
for i in range(len(df)):
	chrom = df.iloc[i].chr
	cary = df.iloc[i].caryotype
	gclust = df.iloc[i].gmm_clust
	sub_dfGBy =  dfGpByClust[(dfGpByClust.chr==chrom)&(dfGpByClust.gmm_clust==gclust)]
	for l in range(len(sub_dfGBy)):
	#df[i:i+1]['final_cluster'] = dfGpByClust[(dfGpByClust.chr==df.iloc[i].chr)&(dfGpByClust.gmm_clust==df.iloc[i].gmm_clust)].K_chosen
		if sub_dfGBy.iloc[l].caryotype == cary :
			titi[i]= sub_dfGBy.iloc[l].K_chosen
			continue

df["final_cluster"] = titi

# Answers
# moyenne ponderee par le nombre de mutations dans cluster

a1 = max(dfGpByClust.groupby(dfGpByClust.K_chosen).apply(lambda x: np.average(x.Pclonal_chosen, weights=x.mut_count)))
b1 = len(dfGpByClust.groupby(dfGpByClust.K_chosen).apply(lambda x: np.average(x.Pclonal_chosen, weights=x.mut_count)))
c1 = pd.DataFrame({
'mut_count_sum' : dfGpByClust.groupby(dfGpByClust.K_chosen)[("mut_count")].sum(),
'pClonal_mean' : dfGpByClust.groupby(dfGpByClust.K_chosen).apply(lambda x: np.average(x.Pclonal_chosen, weights=x.mut_count))})

# write files
#1A
with open(argv[2], 'w') as the_file:
	the_file.write("%f"%a1)
#1B
with open(argv[3], 'w') as the_file:
	the_file.write("%d"%b1)
# 1C
c1.to_csv(argv[4],sep='\t',header=False,index=True)
# 2A
df["final_cluster"].to_csv(argv[5],sep='\t',header=False,index=False)

