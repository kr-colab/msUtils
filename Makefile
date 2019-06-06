CC = gcc
CFLAGS = -O3 -Wall -lm 

msMask:	msMask.c 
	$(CC) msMask.c  -o msMask $(CFLAGS)

msMaskAllRows:	msMaskAllRows.c 
	$(CC) msMaskAllRows.c  -o msMaskAllRows $(CFLAGS)
	
maskedStats:	maskedStats.c msGeneralStats.c
	$(CC) maskedStats.c msGeneralStats.c -o maskedStats $(CFLAGS)
	
niceStatsDan:	niceStatsDan.c msGeneralStats.c
		$(CC) niceStatsDan.c msGeneralStats.c -g -o niceStatsDan $(CFLAGS)

niceStatsSFSRegular:	niceStatsSFSRegular.c msGeneralStats.c
		$(CC) niceStatsSFSRegular.c msGeneralStats.c -g -o niceStatsSFSRegular $(CFLAGS)

niceStats:	niceStats.c msGeneralStats.c
		$(CC) niceStats.c msGeneralStats.c -o niceStats $(CFLAGS)
		
maskedStatsSubpop:	maskedStatsSubpop.c msGeneralStats.c
	$(CC) maskedStatsSubpop.c msGeneralStats.c -o maskedStatsSubpop $(CFLAGS)

twoPopnNiceStats:	twoPopnNiceStats.c msGeneralStats.c
	$(CC) twoPopnNiceStats.c msGeneralStats.c -o twoPopnNiceStats $(CFLAGS)
	
twoPopnStats_forML:	twoPopnStats_forML.c msGeneralStats.c
	$(CC) twoPopnStats_forML.c msGeneralStats.c -o twoPopnStats_forML $(CFLAGS) 	

onePopnStats_forGhostIntroML:	onePopnStats_forGhostIntroML.c msGeneralStats.c
	$(CC) onePopnStats_forGhostIntroML.c msGeneralStats.c -o onePopnStats_forGhostIntroML $(CFLAGS) 	

threePopnStats:	threePopnStats.c msGeneralStats.c
		$(CC) threePopnStats.c msGeneralStats.c -o threePopnStats $(CFLAGS)
msParams: msParams.c
	$(CC) msParams.c ../coalLib/ranlibComplete.c ../pgLib/bedFile.c -o msParams $(CFLAGS)

msParamsSubpop: msParamsSubpop.c
	$(CC) msParamsSubpop.c ../coalLib/ranlibComplete.c ../pgLib/bedFile.c -o msParamsSubpop $(CFLAGS)

msParamsSubpopNoAd: msParamsSubpopNoAd.c
	$(CC) msParamsSubpopNoAd.c ../coalLib/ranlibComplete.c ../pgLib/bedFile.c -o msParamsSubpopNoAd $(CFLAGS)
	
msParamsSubpopTrans: msParamsSubpopTrans.c
	$(CC) msParamsSubpopTrans.c ../coalLib/ranlibComplete.c ../pgLib/bedFile.c -o msParamsSubpopTrans $(CFLAGS)

msParamsTest: msParamsTest.c
	$(CC) msParamsTest.c ../coalLib/ranlibComplete.c ../pgLib/bedFile.c -o msParamsTest $(CFLAGS)

pairDist:	pairDist.c msGeneralStats.c
	$(CC) pairDist.c msGeneralStats.c -o pairDist $(CFLAGS)
pairwiseDists:	pairwiseDists.c msGeneralStats.c
	$(CC) pairwiseDists.c msGeneralStats.c -o pairwiseDists $(CFLAGS)

pairwiseIBSTracts:	pairwiseIBSTracts.c msGeneralStats.c
	$(CC) pairwiseIBSTracts.c msGeneralStats.c -o pairwiseIBSTracts $(CFLAGS)
	
msHKA: msHKA.c msGeneralStats.c
	$(CC) msHKA.c msGeneralStats.c -o msHKA $(CFLAGS)
	
ms2TwoSite: ms2TwoSite.c msGeneralStats.c
	$(CC) ms2TwoSite.c msGeneralStats.c -o ms2TwoSite $(CFLAGS)

ms2RSquareByDistance: ms2RSquareByDistance.c msGeneralStats.c
	$(CC) ms2RSquareByDistance.c msGeneralStats.c -o ms2RSquareByDistance $(CFLAGS)

ms2TwoSite2Popn: ms2TwoSite2Popn.c msGeneralStats.c
	$(CC) ms2TwoSite2Popn.c msGeneralStats.c -o ms2TwoSite2Popn $(CFLAGS)

ms2SFS2D: ms2SFS2D.c msGeneralStats.c
	$(CC) ms2SFS2D.c msGeneralStats.c -o ms2SFS2D $(CFLAGS)
ms2SFS2D_discoal: ms2SFS2D_discoal.c msGeneralStats.c
	$(CC) ms2SFS2D_discoal.c msGeneralStats.c -o ms2SFS2D_discoal $(CFLAGS)

ms2SFSVector: ms2SFSVector.c msGeneralStats.c
	$(CC) ms2SFSVector.c msGeneralStats.c -o ms2SFSVector $(CFLAGS)

ms2SFSVectorWindow: ms2SFSVectorWindow.c msGeneralStats.c
	$(CC) ms2SFSVectorWindow.c msGeneralStats.c -o ms2SFSVectorWindow $(CFLAGS)
	
discoal_mig2hmm: discoal_mig2hmm.c msGeneralStats.c
	$(CC) discoal_mig2hmm.c msGeneralStats.c -o discoal_mig2hmm $(CFLAGS)

slideFST: slideFST.c msGeneralStats.c
	$(CC) slideFST.c msGeneralStats.c -o slideFST $(CFLAGS)

niceStatsNoOmega:	niceStatsNoOmega.c msGeneralStats.c
		$(CC) niceStatsNoOmega.c msGeneralStats.c -g -o niceStatsNoOmega $(CFLAGS)	

niceStatsShanku:	niceStatsShanku.c msGeneralStats.c
		$(CC) niceStatsShanku.c msGeneralStats.c -o niceStatsShanku $(CFLAGS)

niceStatsAchazSystem:	niceStatsAchazSystem.c msGeneralStats.c
		$(CC) niceStatsAchazSystem.c msGeneralStats.c -O3 -o niceStatsAchazSystem $(CFLAGS)
		
niceStatsDiploid:	niceStatsDiploid.c msGeneralStats.c msDiploidStats.c
		$(CC) niceStatsDiploid.c msGeneralStats.c msDiploidStats.c -o niceStatsDiploid $(CFLAGS)
niceStats4Gamete:	niceStats4Gamete.c msGeneralStats.c 
			$(CC) niceStats4Gamete.c msGeneralStats.c -o niceStats4Gamete $(CFLAGS)	

