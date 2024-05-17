library(BrainTransR)


tmp=get_geneSetList(type='GO',parameter = 'BP')
tmp=get_geneSetList(type='GO',parameter = 'CC')
tmp=get_geneSetList(type='GO',parameter = 'MF')
tmp=get_geneSetList(type='KEGG')
tmp=get_geneSetList(type='Reactome')
tmp=get_geneSetList(type='WikiPathways')
tmp=get_geneSetList(type='SynGO')
tmp=get_geneSetList(type='CellType', parameter = 'Seidlitz2020')
tmp=get_geneSetList(type='CellType', parameter = 'Lake2018')
tmp=get_geneSetList(type='CellType', parameter = 'Martins2021')



