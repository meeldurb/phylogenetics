##install.packages('circlize',  repos = "http://cran.rstudio.com/" )  
#install.packages('data.table', repos = "http://cran.rstudio.com/")

library(circlize)
library(data.table)




karyo <- fread(colClasses = c("character", "numeric",
                              "numeric", "character", "character"),
               "ssa01  0          45670927   1p    gneg
               ssa01  45670927   102921759  1qa   gpos50
               ssa01  102921759  158994226  1qb   gneg
               ssa02  0          37194884   2p    gpos50
               ssa02  37194884   72887026   2q    gneg
               ssa03  0          53325515   3p    gpos50
               ssa03  53325515   92451775   3q    gneg
               ssa04  0          24769547   4p    gpos50
               ssa04  24769547   82368206   4q    gneg
               ssa05  0          44366638   5p    gpos50
               ssa05  44366638   80452451   5q    gneg
               ssa06  0          44681795   6p    gpos50
               ssa06  44681795   87037992   6q    gneg
               ssa07  0          29822415   7p    gpos50
               ssa07  29822415   58755209   7q    gneg
               ssa08  -15000000  0          rDNA  gpos50
               ssa08  0          26442203   8q    gneg
               ssa09  0          48560009   9qa   gpos50
               ssa09  48560009   92800824   9qb   gneg
               ssa09  92800824   141678924  9qc   gpos50
               ssa10  0          50712303   10qa  gneg
               ssa10  50712303   116119932  10qb  gpos50
               ssa11  0          44252459   11qa  gneg
               ssa11  44252459   93835980   11qb  gpos50
               ssa12  0          31366442   12qa  gneg
               ssa12  31366442   91841354   12qb  gpos50
               ssa13  0          46825690   13qa  gneg
               ssa13  46825690   107755956  13qb  gpos50
               ssa14  0          45561258   14qa  gneg
               ssa14  45561258   93876043   14qb  gpos50
               ssa15  0          53413608   15qa  gneg
               ssa15  53413608   103934464  15qb  gpos50
               ssa16  0          56974611   16qa  gneg
               ssa16  56974611   87752276   16qb  gpos50
               ssa17  0          31608156   17qa  gneg
               ssa17  31608156   57636407   17qb  gpos50
               ssa18  0          35154609   18qa  gneg
               ssa18  35154609   70656513   18qb  gpos50
               ssa19  0          38868197   19qa  gneg
               ssa19  38868197   82958981   19qb  gpos50
               ssa20  0          40962200   20qa  gneg
               ssa20  40962200   86771154   20qb  gpos50
               ssa21  0  58011450  21  gneg
               ssa22  0  63437620  22  gneg
               ssa23  0  49935705  23  gneg
               ssa24  0  48640078  24  gneg
               ssa25  0  51465665  25  gneg
               ssa26  0  47887414  26  gneg
               ssa27  0  43932925  27  gneg
               ssa28  0  39590257  28  gneg
               ssa29  0  42472244  29  gneg")

cells <- fread("ssa01	1614959	1669305	ssa09	29281275	29333044	cigene_21_1
ssa01	7016802	7020787	ssa09	4823610	4828206	cigene_21_1
               ssa01	9081979	9093865	ssa09	154581	166196	cigene_21_1
               ssa01	9712460	9739135	ssa09	7445939	7465960	cigene_21_1
               ssa01	12071020	12408277	ssa09	9274905	9632122	cigene_21_1
               ssa01	14572869	42474127	ssa09	11608030	41766237	cigene_21_1
               ssa01	46210305	47454891	ssa18	3719997	4493435	cigene_21_1
               ssa01	50060796	50239648	ssa18	7424831	7583455	cigene_21_1
               ssa01	51210061	51300453	ssa18	8545868	8626755	cigene_21_1
               ssa01	51897468	56353439	ssa18	9396892	13754255	cigene_21_1
               ssa01	58273089	71773335	ssa18	15568806	28115019	cigene_21_1
               ssa01	75701612	84617654	ssa28	15206662	21908859	cigene_21_1
               ssa01	88946361	89904339	ssa18	29057097	30411506	cigene_21_1
               ssa01	91244470	98677754	ssa28	25483092	32239461	cigene_21_1
               ssa01	109309626	121051890	ssa13	86130521	97835264	cigene_21_1
               ssa01	122036199	121252179	ssa11	71827408	72738102	cigene_21_1
               ssa01	122475995	122615318	ssa11	71193623	71354981	cigene_21_1
               ssa01	125307412	146930885	ssa11	46241584	68505601	cigene_21_1
               ssa01	153271788	153799816	ssa13	101270545	101940810	cigene_21_1
               ssa02	11684996	12717949	ssa05	67371143	69170152	cigene_21_2
               ssa02	16546883	16555231	ssa05	64129221	64141747	cigene_21_2
               ssa02	18642703	36708012	ssa05	45090772	62193486	cigene_21_2
               ssa02	37941652	47370186	ssa12	21170522	29959482	cigene_21_2
               ssa02	47620408	47680373	ssa12	20878203	20940381	cigene_21_2
               ssa02	48081545	48089303	ssa12	20427894	20435961	cigene_21_2
               ssa02	52121768	52128448	ssa12	16669348	16754965	cigene_21_2
               ssa03	5220495	5278171	ssa14	7188182	7300413	cigene_21_3
               ssa03	6293360	6305698	ssa14	5284513	5297392	cigene_21_3
               ssa03	7505213	7516846	ssa14	6444195	6457191	cigene_21_3
               ssa03	8076589	8278973	ssa14	7785768	7996369	cigene_21_3
               ssa03	9529043	9594449	ssa14	8736981	8786214	cigene_21_3
               ssa03	9989820	9992111	ssa14	9188143	9190512	cigene_21_3
               ssa03	10036701	36406504	ssa14	9238821	35812404	cigene_21_3
               ssa03	37416988	41498999	ssa14	37354027	41614602	cigene_21_3
               ssa03	43523350	44059872	ssa14	43140266	43637637	cigene_21_3
               ssa03	49493335	49516186	ssa23	1125100	1154182	cigene_21_3
               ssa03	49934274	49989103	ssa23	1473274	1512798	cigene_21_3
               ssa03	51705175	52164612	ssa23	2622321	3091316	cigene_21_3
               ssa03	52657938	52671218	ssa23	156538	161076	cigene_21_3
               ssa03	53095574	53156237	ssa23	553450	643340	cigene_21_3
               ssa03	53401956	60247561	ssa06	38055328	44367386	cigene_21_3
               ssa03	60804564	60625652	ssa06	37334436	37542208	cigene_21_3
               ssa03	63510807	63516671	ssa06	34703519	34710646	cigene_21_3
               ssa03	65313004	78747949	ssa06	21290612	33904873	cigene_21_3
               ssa04	15892514	15916612	ssa08	577512	598283	cigene_21_4
               ssa04	16094947	16102388	ssa08	859734	867857	cigene_21_4
               ssa04	16691555	24323800	ssa08	1125618	7557097	cigene_21_4
               ssa04	31343756	39113764	ssa11	72904946	79771408	cigene_21_4
               ssa04	39475843	77327371	ssa13	46934906	85580013	cigene_21_4
               ssa05	7238640	7255734	ssa09	85035044	85056629	cigene_21_5
               ssa05	9764469	9770104	ssa09	82842165	82848632	cigene_21_5
               ssa05	11381815	11396163	ssa09	79954254	79962499	cigene_21_5
               ssa05	11734902	11736744	ssa09	80647653	80650012	cigene_21_5
               ssa05	12406685	12412632	ssa09	79177037	79182329	cigene_21_5
               ssa05	12543696	12549047	ssa09	83219885	83224777	cigene_21_5
               ssa05	20847234	20926816	ssa09	49355287	49508142	cigene_21_5
               ssa05	22644558	40694462	ssa09	51558992	71436753	cigene_21_5
               ssa06	45455711	63140141	ssa15	22395252	38332623	cigene_21_6
               ssa06	64749858	78264252	ssa15	42108401	53651273	cigene_21_6
               ssa06	80383575	80403209	ssa15	19148552	19150469	cigene_21_6
               ssa06	81626927	81631894	ssa15	16988920	16991022	cigene_21_6
               ssa07	6963739	7160574	ssa18	42582275	42797976	cigene_21_7
               ssa07	7993826	8059425	ssa18	44228520	44314804	cigene_21_7
               ssa07	9017657	9021090	ssa18	45079407	45084561	cigene_21_7
               ssa07	9727784	9739801	ssa18	45831710	45846205	cigene_21_7
               ssa07	11707836	16557135	ssa18	47611914	52852095	cigene_21_7
               ssa07	18853169	25142490	ssa18	55094003	64180758	cigene_21_7
               ssa07	30811148	45673344	ssa17	31983813	47865264	cigene_21_7
               ssa09	92442544	129272924	ssa20	404265	73606032	cigene_21_8
               ssa10	8167871	8282101	ssa16	52932419	52991335	cigene_21_9
               ssa10	8766821	8793303	ssa16	52434215	52458565	cigene_21_9
               ssa10	38218743	10018175	ssa16	12248792	49674543	cigene_21_9
               ssa10	43022698	44205862	ssa23	7755007	8483599	cigene_21_9
               ssa10	46906723	47424432	ssa23	4509729	4933463	cigene_21_9
               ssa10	48752846	48810277	ssa23	3681017	773826	cigene_21_9
               ssa10	65752888	72629719	ssa23	34814178	41802097	cigene_21_9
               ssa10	73929991	101902826	ssa16	10071611	37290098	cigene_21_9
               ssa10	102512162	102528434	ssa16	9486446	9533380	cigene_21_9
               ssa10	102829653	102841157	ssa16	9114869	9136730	cigene_21_9
               ssa10	105335381	106742980	ssa16	4811302	6310072	cigene_21_9
               ssa10	106901912	106929882	ssa16	5193702	5241920	cigene_21_9
               ssa10	107625873	108769684	ssa16	3212317	4467957	cigene_21_9
               ssa10	112648744	112848944	ssa23	44161388	44432184	cigene_21_9
               ssa11	6947124	9464323	ssa26	6838726	9499944	cigene_21_10
               ssa11	9756254	14778105	ssa26	9801433	10287714	cigene_21_10
               ssa11	14818665	21376010	ssa26	15198949	22351654	cigene_21_10
               ssa11	24165857	24180894	ssa26	30792247	30805890	cigene_21_10
               ssa11	24482078	27092054	ssa26	26544906	29238954	cigene_21_10
               ssa12	34790044	34794647	ssa22	387405	399983	cigene_21_11
               ssa12	82171329	34980877	ssa22	5347453	55995590	cigene_21_11
               ssa12	88947419	88950644	ssa22	58061857	58065885	cigene_21_11
               ssa13	4696542	5407469	ssa15	97165749	98128179	cigene_21_12
               ssa13	6817935	6868364	ssa15	96354751	96407117	cigene_21_12
               ssa13	7859059	9945959	ssa15	71333556	73138997	cigene_21_12
               ssa13	16111515	22078752	ssa15	56277269	62732856	cigene_21_12
               ssa13	22341332	25690806	ssa15	75842685	79443004	cigene_21_12
               ssa13	26414630	27800548	ssa15	64133753	65632535	cigene_21_12
               ssa13	28988629	31445267	ssa15	66117700	68360923	cigene_21_12
               ssa13	31910188	42698339	ssa15	80352580	93307813	cigene_21_12
               ssa13	45856726	45913743	ssa15	69588302	69642314	cigene_21_12
               ssa14	52839030	54384200	ssa27	1443466	3184433	cigene_21_13
               ssa14	56821525	57152716	ssa27	7724273	8391065	cigene_21_13
               ssa14	57493930	71624769	ssa27	8055511	22518845	cigene_21_13
               ssa14	72082471	75739508	ssa27	32294210	35932242	cigene_21_13
               ssa14	78558507	83963615	ssa27	24094384	29533768	cigene_21_13
               ssa14	86158407	86171197	ssa27	36904661	36920410	cigene_21_13
               ssa15	1752013	1770993	ssa24	435977	525434	cigene_21_14
               ssa15	4083941	4089893	ssa24	43151231	43156913	cigene_21_14
               ssa16	59189479	59380167	ssa17	1657748	1890766	cigene_21_15
               ssa16	60965863	61091594	ssa17	3574607	3696442	cigene_21_15
               ssa16	62432999	62738867	ssa17	5113928	5559321	cigene_21_15
               ssa16	63946386	72264032	ssa17	6697997	14684669	cigene_21_15
               ssa18	1421949	1569141	ssa19	65140287	65272500	cigene_21_16
               ssa18	2096386	2237923	ssa19	64607220	64731693	cigene_21_16
               ssa18	2597795	2792841	ssa19	64206783	64148327	cigene_21_16
               ssa19	8211366	8224304	ssa29	27586448	27594245	cigene_21_17
               ssa19	10755365	13926936	ssa29	31423272	35634434	cigene_21_17
               ssa19	18740368	18780885	ssa29	2218060	2292146	cigene_21_17
               ssa19	19649693	19757868	ssa29	3641401	3752843	cigene_21_17
               ssa19	20377463	20799530	ssa29	4522106	4635725	cigene_21_17
               ssa19	21082712	21335242	ssa29	5061216	5284027	cigene_21_17
               ssa19	22000576	22071595	ssa29	5641744	5688443	cigene_21_17
               ssa19	23505136	23588019	ssa29	6683526	6789047	cigene_21_17
               ssa19	23813680	23814997	ssa29	7066961	7069972	cigene_21_17
               ssa19	24953824	24970036	ssa29	8424176	8435402	cigene_21_17
               ssa19	26036962	26309592	ssa29	9318270	9598400	cigene_21_17
               ssa19	26762459	31928257	ssa29	10053635	14057287	cigene_21_17
               ssa19	38596874	38853525	ssa28	34744018	35014266	cigene_21_17
               ssa19	39628318	46358256	ssa28	6210175	12239122	cigene_21_17
               ssa19	47372529	47747106	ssa28	4605505	5096360	cigene_21_17
               ssa19	48414218	50045701	ssa28	2231796	3867555	cigene_21_17
               ssa19	50448845	50694195	ssa28	1424461	1725818	cigene_21_17
               ssa19	50923185	50992403	ssa28	1033117	1097148	cigene_21_17
               ssa19	51131001	51160567	ssa28	746790	770334	cigene_21_17
               ssa19	51325981	51337975	ssa28	543837	571363	cigene_21_17
               ssa19	66573521	75061073	ssa29	19305856	27238917	cigene_21_17
               ssa20	1103917	2131633	ssa24	3507761	4623923	cigene_21_18
               ssa20	6164522	13874219	ssa24	31190819	39448466	cigene_21_18
               ssa20	16769388	40764837	ssa24	26806440	4985055	cigene_21_18
               ssa21	7434216	16880593	ssa25	5293552	14686352	cigene_21_19
               ssa21	16872721	16880593	ssa25	14679865	14686352	cigene_21_19
               ssa21	19117340	50737453	ssa25	15339084	44368707	cigene_21_19
               ")
names(cells) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "col")
cells <- data.frame(cells)

regions_annot <- fread("ssa01  0          10000000   ssa01  region_telomeric
                       ssa01  10000000   148994226  ssa01  region_normal
                       ssa01  148994226  158994226  ssa01  region_telomeric
                       ssa02  0          12000000   ssa02  region_collapsed
                       ssa02  12000001   51206500   ssa02  region_highid
                       ssa02  51206500   72887026   ssa02  region_collapsed
                       ssa03  0          10000000   ssa03  region_telomeric
                       ssa03  10000000   50048600   ssa03  region_normal
                       ssa03  50048600   69337097   ssa03  region_highid
                       ssa03  69337097   92451775   ssa03  region_collapsed
                       ssa04  0          11000000   ssa04  region_collapsed
                       ssa04  11000001   25621004   ssa04  region_highid
                       ssa04  25621004   72368206   ssa04  region_normal
                       ssa04  72368206   82368206   ssa04  region_telomeric
                       ssa05  0          8000000    ssa05  region_highid
                       ssa05  8000000    43452452   ssa05  region_normal
                       ssa05  43452452   69452452   ssa05  region_highid
                       ssa05  69452452   80452451   ssa05  region_collapsed
                       ssa06  0          27497429   ssa06  region_collapsed
                       ssa06  27497429   44866161   ssa06  region_highid
                       ssa06  44866161   77037992   ssa06  region_normal
                       ssa06  77037992   87037992   ssa06  region_telomeric
                       ssa07  0          10000000   ssa07  region_telomeric
                       ssa07  10000000   31099723   ssa07  region_normal
                       ssa07  31099723   47099723   ssa07  region_highid
                       ssa07  47099723   58755209   ssa07  region_collapsed
                       ssa08  0          16055992   ssa08  region_highid
                       ssa08  16055992   26442203   ssa08  region_collapsed
                       ssa09  0          84697819   ssa09  region_normal
                       ssa09  84697819   141678924  ssa09  region_highid
                       ssa10  0          106119932  ssa10  region_normal
                       ssa10  106119932  116119932  ssa10  region_telomeric
                       ssa11  0          6847202    ssa11  region_normal
                       ssa11  6847202    28804710   ssa11  region_highid
                       ssa11  28804710   44804709   ssa11  region_collapsed
                       ssa11  44804709   83835980   ssa11  region_normal
                       ssa11  83835980   93835980   ssa11  region_telomeric
                       ssa12  0          16357447   ssa12  region_collapsed
                       ssa12  16357447   30357447   ssa12  region_highid
                       ssa12  30357447   81841354   ssa12  region_normal
                       ssa12  81841354   91841354   ssa12  region_telomeric
                       ssa13  0          97755956   ssa13  region_normal
                       ssa13  97755956   107755956  ssa13  region_telomeric
                       ssa14  0          83876043   ssa14  region_normal
                       ssa14  83876043   93876043   ssa14  region_telomeric
                       ssa15  0          93934464   ssa15  region_normal
                       ssa15  93934464   103934464  ssa15  region_telomeric
                       ssa16  0          58868232   ssa16  region_normal
                       ssa16  58868232   73868231   ssa16  region_highid
                       ssa16  73868231   87752276   ssa16  region_collapsed
                       ssa17  0          16168292   ssa17  region_highid
                       ssa17  16168292   24168292   ssa17  region_collapsed
                       ssa17  24168293   48965773   ssa17  region_highid
                       ssa17  48965773   57636407   ssa17  region_collapsed
                       ssa18  0          60656513   ssa18  region_normal
                       ssa18  60656513   70656513   ssa18  region_telomeric
                       ssa19  0          72958981   ssa19  region_normal
                       ssa19  72958981   82958981   ssa19  region_telomeric
                       ssa20  0          40968350   ssa20  region_normal
                       ssa20  40968350   86771154   ssa20  region_highid
                       ssa21  0          48011450   ssa21  region_normal
                       ssa21  48011450   58011450   ssa21  region_telomeric
                       ssa22  0          53437620   ssa22  region_normal
                       ssa22  53437620   63437620   ssa22  region_telomeric
                       ssa23  0          39935705   ssa23  region_normal
                       ssa23  39935705   49935705   ssa23  region_telomeric
                       ssa24  0          38640078   ssa24  region_normal
                       ssa24  38640078   48640078   ssa24  region_telomeric
                       ssa25  0          41465665   ssa25  region_normal
                       ssa25  41465665   51465665   ssa25  region_telomeric
                       ssa26  0          6517046    ssa26  region_normal
                       ssa26  6517046    29138210   ssa26  region_highid
                       ssa26  29138210   47887414   ssa26  region_collapsed
                       ssa27  0          33932925   ssa27  region_normal
                       ssa27  33932925   43932925   ssa27  region_telomeric
                       ssa28  0          29590257   ssa28  region_normal
                       ssa28  29590257   39590257   ssa28  region_telomeric
                       ssa29  0          32472244   ssa29  region_normal
                       ssa29  32472244   42472244   ssa29  region_telomeric")
names(regions_annot) <- c("chr", "start", "end", "value", "col")

getCol <- function(r, g, b, alpha=0.8) {
  rgb(r/255, g/255, b/255, alpha=alpha)
}

getColName <- function(col.in, alpha=0.8) {
  co = as.numeric(col2rgb(col.in))
  rgb(red=co[1]/255, green=co[2]/255, blue=co[3]/255, alpha=alpha)
}

getColFromList <- function(color) {
  if (!is.function(get(color))) { return(getColName(color)) } # "grey", "black"
  else { return(getColName(color)) } #cigene_21_1"
}

alpha <- 0.8
cigene_21_1 <- getCol(221,119,136,alpha)
cigene_21_12 <- getCol(170,68,136,alpha)
cigene_21_13 <- getCol(204,153,187,alpha)
cigene_21_18 <- getCol(17,68,119,alpha)
cigene_21_3 <- getCol(68,119,170,alpha)
cigene_21_6 <- getCol(119,170,221,alpha)
cigene_21_5 <- getCol(17,119,119,alpha)
cigene_21_16 <- getCol(170,68,85,alpha)
cigene_21_8 <- getCol(119,204,204,alpha)
cigene_21_21 <- getCol(17,119,68,alpha)
cigene_21_9 <- getCol(68,170,119,alpha)
cigene_21_10 <- getCol(136,204,170,alpha)
cigene_21_19 <- getCol(119,119,17,alpha)
cigene_21_14 <- getCol(170,170,68,alpha)
cigene_21_17 <- getCol(221,221,119,alpha)
cigene_21_2 <- getCol(119,68,17,alpha)
cigene_21_15 <- getCol(170,119,68,alpha)
cigene_21_4 <- getCol(221,170,119,alpha)
cigene_21_11 <- getCol(119,17,34,alpha)
cigene_21_20 <- getCol(68,170,170,alpha)
cigene_21_7 <- getCol(119,17,85,alpha)

alpha <- 1
region_telomeric <- getCol(255,255,191,alpha)
region_normal <- getCol(102,194,165,alpha)
region_collapsed <- getCol(213,62,79,alpha)
region_highid <- getCol(244,109,67,alpha)


regions_annot.bed <- data.frame(regions_annot)[1:3]
regions_annot.bed$value1 <- rep(1, nrow(regions_annot.bed))
regions_annot.bed$value2 <- rep(0, nrow(regions_annot.bed))
regions_annot.bed$col <- as.character(unlist(mget(regions_annot$col, ifnotfound = list(function(x) getColFromList(x)))))

cells$col <- as.character(unlist(mget(cells$col, ifnotfound = list(function(x) getColFromList(x)))))

par(mar = c(1, 1, 1, 1))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(data.frame(karyo), track.height=0.05, ideogram.height = 0.03, 
                              chromosome.index = sort(unique(karyo$V1)))
circos.genomicTrackPlotRegion(regions_annot.bed, track.height=0.02, 
                              panel.fun = function(region,  value,  ...) {
  circos.genomicRect(region, value, col = value$col, border = NA, ...)
})

Ss4R.time.chr <- read.table(file = "20170919-Ss4r_time&chropos_fortrack.csv", 
                            sep = ";", stringsAsFactors = F, header = T)
Ss4R.time.chr <- na.omit(Ss4R.time.chr)


circos.genomicTrackPlotRegion(data.frame(Ss4R.time.chr[,1:4]), track.height=0.2, 
                              panel.fun = function(region, value,...){
  circos.genomicPoints(region, value, type="l", col="dark red", cex = 0.2, border=NA, ylim=c(0,1500),...)
})
circos.genomicLink(cells[,1:3], cells[,4:6], col = cells$col, border = NA)

