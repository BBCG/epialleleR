test_generateCytosineReport <- function () {
  amplicon.bam    <- system.file("extdata", "amplicon010meth.bam", package="epialleleR")
  amplicon.filter <- generateCytosineReport(amplicon.bam, verbose=TRUE)
  amplicon.nofilter <- generateCytosineReport(amplicon.bam, filter.reads=FALSE, verbose=TRUE)
  amplicon.nothresh <- generateCytosineReport(amplicon.bam, threshold.reads=FALSE, verbose=TRUE)
  amplicon.donothing <- generateCytosineReport(amplicon.bam, filter.reads=FALSE, threshold.reads=FALSE, verbose=TRUE)
  
  RUnit::checkEquals(
    amplicon.filter[, c(sum(meth), sum(unmeth))],
    c(632, 6430)
  )
  RUnit::checkEquals(
    amplicon.nofilter[, c(sum(meth), sum(unmeth))],
    c(636, 6445)
  )
  RUnit::checkEquals(
    amplicon.nothresh[, c(sum(meth), sum(unmeth))],
    c(679, 6383)
  )
  RUnit::checkEquals(
    amplicon.donothing[, c(sum(meth), sum(unmeth))],
    c(683, 6398)
  )
  
  capture.bam <- system.file("extdata", "capture.bam", package="epialleleR")
  cg.report   <- generateCytosineReport(capture.bam, verbose=TRUE)
  cx.report   <- generateCytosineReport(capture.bam, threshold.reads=FALSE,
                                        report.context="CX", verbose=FALSE)
  cx.nofilt   <- generateCytosineReport(capture.bam, threshold.reads=FALSE,
                                        filter.reads=FALSE, report.context="CX")
  
  RUnit::checkEquals(
    nrow(cx.report[duplicated(paste(rname,pos))]),
    0
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.report$strand)[c("+", "-")]),
    c(41522, 41233)
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.nofilt$strand)[c("+", "-")]),
    c(48517, 48669)
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.report$context)[c("CHH", "CHG", "CG")]),
    c(47319, 20269, 15167)
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.nofilt$context)[c("CHH", "CHG", "CG")]),
    c(58292, 23486, 15408)
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.report[strand=="+"]$context)[c("CHH", "CHG", "CG")]),
    c(23515, 10222, 7785)
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.nofilt[strand=="+"]$context)[c("CHH", "CHG", "CG")]),
    c(28762, 11853, 7902)
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.report[strand=="-"]$context)[c("CHH", "CHG", "CG")]),
    c(23804, 10047, 7382)
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.nofilt[strand=="-"]$context)[c("CHH", "CHG", "CG")]),
    c(29530, 11633, 7506)
  )
  
  RUnit::checkEquals(
    dim(cg.report),
    c(15167,6)
  )
  
  RUnit::checkEquals(
    dim(cx.report),
    c(82755,6)
  )
  
  RUnit::checkEquals(
    dim(cx.nofilt),
    c(97186,6)
  )
  
  RUnit::checkEquals(
    sum(cg.report$meth),
    4974
  )

  RUnit::checkEquals(
    sum(cg.report$unmeth),
    14859
  )
  
  RUnit::checkEquals(
    sum(cx.report$meth),
    5732
  )
  
  RUnit::checkEquals(
    sum(cx.report$unmeth),
    103454
  )
  
  RUnit::checkEquals(
    sum(cx.nofilt$meth),
    6051
  )
  
  RUnit::checkEquals(
    sum(cx.nofilt$unmeth),
    125903
  )
  
  # some extended consistency checks
  RUnit::checkEquals(
    cx.nofilt[context=="CG", sum(meth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(222, 242, 128, 91, 167, 172, 101, 77, 85, 18, 81, 64, 159, 240, 116, 105, 129, 140, 16, 39, 107, 81,
      161, 62, 59, 31, 140, 104, 73, 37, 181, 103, 406, 457, 13, 4, 63, 90, 253, 438, 91, 56, 15, 22, 106, 91)
  )
  RUnit::checkEquals(
    cx.nofilt[context=="CG", sum(unmeth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(556, 713, 295, 316, 446, 679, 162, 115, 109, 82, 404, 289, 256, 336, 341, 326, 326, 102, 328, 207, 446, 609, 328,
      180, 148, 163, 243, 267, 283, 198, 535, 482, 1035, 1064, 97, 80, 177, 270, 447, 394, 65, 37, 92, 157, 197, 201)
  )
  RUnit::checkEquals(
    cx.nofilt[context=="CHG", sum(meth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(4, 6, 2, 1, 5, 7, 2, 1, 1, 2, 4, 1, 1, 2, 3, 4, 3, 4, 4, 1, 2, 4, 1, 1, 5, 2, 3, 0, 1, 2, 2, 5, 11, 5, 1, 0, 1, 1, 8, 3, 1, 2, 0, 1, 1, 2)
  )
  RUnit::checkEquals(
    cx.nofilt[context=="CHG", sum(unmeth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(1149, 1258, 767, 693, 937, 1204, 349, 312, 291, 226, 754, 524, 616, 819, 683, 849, 825, 472, 393, 349, 821, 955, 751,
      462, 301, 270, 582, 546, 574, 384, 1267, 1063, 2488, 2507, 110, 67, 448, 496, 883, 1287, 285, 191, 150, 191, 595, 550)
  )
  RUnit::checkEquals(
    cx.nofilt[context=="CHH", sum(meth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(9, 14, 6, 6, 14, 12, 1, 4, 3, 4, 8, 6, 3, 10, 9, 8, 6, 6, 2, 4, 4, 10, 6, 3, 3, 5, 6, 0, 4, 6, 10, 9, 26, 17, 1, 0, 3, 7, 11, 8, 2, 3, 3, 0, 5, 5)
  )
  RUnit::checkEquals(
    cx.nofilt[context=="CHH", sum(unmeth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(2921, 3293, 1615, 1589, 2415, 3348, 735, 983, 720, 674, 1824, 1306, 1495, 1886, 2008, 2295, 1925, 1366, 950, 850, 2199, 2487, 1732,
      1026, 846, 829, 1567, 1212, 1249, 956, 2675, 2522, 6963, 6300, 304, 204, 968, 1350, 2144, 2775, 620, 513, 344, 371, 1634, 1638)
  )
  
  
  generateCytosineReport(capture.bam, report.file=tempfile())
  
  cx.trim <- generateCytosineReport(
    capture.bam, threshold.reads=FALSE, filter.reads=FALSE,
    trim=3, report.context="CX", verbose=FALSE)
  cx.notrim <- generateCytosineReport(
    capture.bam, threshold.reads=FALSE, filter.reads=FALSE,
    trim=0, report.context="CX", verbose=FALSE)
  RUnit::checkTrue(
    ! identical(cx.trim, cx.notrim)
  )
  
  RUnit::checkTrue(
    identical(
      data.table::merge.data.table(
        cx.trim[, .(rname, strand, pos, context)],
        cx.trim[, .(rname, strand, pos, context)]
      ),
      data.table::merge.data.table(
        cx.trim[,   .(rname, strand, pos, context)],
        cx.notrim[, .(rname, strand, pos, context)]
      )
    )
  )
  
  
  cg.quality  <- generateCytosineReport(capture.bam, verbose=TRUE,
                                        min.mapq=30, min.baseq=20)
  cx.quality  <- generateCytosineReport(capture.bam, threshold.reads=FALSE,
                                        filter.reads=FALSE, 
                                        min.mapq=30, min.baseq=20,
                                        report.context="CX", verbose=FALSE)
  cx.filt     <- generateCytosineReport(capture.bam, threshold.reads=FALSE,
                                        filter.reads=TRUE, min.context.sites=5,
                                        max.outofcontext.beta=0.1,
                                        min.mapq=30, min.baseq=20,
                                        report.context="CX", verbose=FALSE)
  
  RUnit::checkEquals(
    dim(cg.quality),
    c(14947,6)
  )
  
  RUnit::checkEquals(
    dim(cx.quality),
    c(96151,6)
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.quality$context)[c("CHH", "CHG", "CG")]),
    c(57687, 23267, 15197)
  )
  
  RUnit::checkEquals(
    sum(cg.quality$meth),
    4830
  )
  
  RUnit::checkEquals(
    sum(cg.quality$unmeth),
    14670
  )
  
  RUnit::checkEquals(
    sum(cx.quality$meth),
    5873
  )
  
  RUnit::checkEquals(
    sum(cx.quality$unmeth),
    124333
  )
  
  # more extended consistency checks
  RUnit::checkEquals(
    cx.quality[context=="CG", sum(meth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(217, 235, 119, 87, 160, 168, 99, 77, 82, 18, 80, 62, 158, 231, 110, 104, 127, 134, 15, 36, 106, 79,
      161, 62, 59, 31, 137, 98, 73, 36, 167, 102, 392, 446, 12, 4, 63, 86, 240, 429, 86, 56, 15, 22, 101, 88)
  )
  RUnit::checkEquals(
    cx.quality[context=="CG", sum(unmeth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(551, 707, 293, 312, 445, 673, 157, 113, 106, 82, 397, 285, 256, 329, 339, 325, 323, 102, 328, 207, 443, 605, 324,
      179, 146, 163, 237, 262, 276, 196, 527, 470, 1027, 1048, 97, 80, 172, 270, 445, 390, 62, 37, 89, 154, 195, 198)
  )
  RUnit::checkEquals(
    cx.quality[context=="CG", sum(as.numeric(pos)), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(81351176333, 73001003519, 34067775901, 32050020375, 45535693217, 54595528147, 31428140625, 22569317363, 24338733458, 10383726558,
      34593702857, 22194464580, 22569211162, 33074156654, 24416818659, 23638910883, 33997597033, 16612624644, 20741550105, 13386587843,
      31190112569, 39171042572, 30154205389, 18624248043, 7078530927, 3948042625, 19450144807, 19693203147, 21553889829, 13527792443,
      24401350957, 23373768915, 41974178009, 35838793619, 3496426056, 3234415920, 7023486782, 12427350607, 11485040910, 11479920727,
      3720934776, 2918776285, 2969393048, 4471621433, 19381406469, 19518705607)
  )
  RUnit::checkEquals(
    cx.quality[context=="CHG", sum(meth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(4, 6, 2, 1, 5, 7, 1, 1, 1, 2, 4, 1, 1, 2, 3, 4, 3, 4, 4, 1, 2, 3, 1, 1, 5, 2, 3, 0, 1, 1, 2, 5, 11, 4, 0, 0, 1, 1, 8, 3, 1, 2, 0, 1, 1, 2)
  )
  RUnit::checkEquals(
    cx.quality[context=="CHG", sum(unmeth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(1131, 1246, 764, 680, 927, 1185, 343, 309, 288, 223, 744, 518, 611, 814, 675, 827, 811, 467, 387, 347, 814, 946, 740,
      458, 298, 266, 579, 536, 569, 377, 1254, 1047, 2459, 2474, 109, 66, 446, 492, 873, 1271, 280, 188, 150, 186, 590, 544)
  )
  RUnit::checkEquals(
    cx.quality[context=="CHG", sum(as.numeric(pos)), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(115471836976, 100376311443, 62233895957, 64390850744, 68919572153, 79885841616, 35996845842, 28726440473, 29993885241, 20239819097,
      54355256540, 35683536419, 33888246497, 48173311519, 42256595388, 50154804462, 59488065292, 33029014705, 26347624033, 22706903579,
      48446890884, 55237929170, 38597655766, 32057668636, 8987612180, 8351801813, 31164951173, 31373834447, 32953932629, 22327010308,
      42917775709, 40377421217, 64872871219, 66761918639, 4002633421, 2505374296, 12497345524, 18843677507, 14376249400, 16016927846,
      6984814573, 5556698053, 4214094025, 4725981470, 41542780004, 37152853172)
  )
  RUnit::checkEquals(
    cx.quality[context=="CHH", sum(meth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(9, 14, 6, 6, 14, 11, 1, 4, 3, 3, 8, 6, 3, 9, 9, 8, 6, 5, 2, 4, 4, 10, 6, 3, 3, 5, 6, 0, 4, 6, 10, 9, 24, 17, 1, 0, 3, 7, 10, 8, 2, 3, 3, 0, 5, 5)
  )
  RUnit::checkEquals(
    cx.quality[context=="CHH", sum(unmeth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(2879, 3239, 1594, 1573, 2402, 3308, 719, 963, 706, 661, 1803, 1287, 1485, 1851, 1982, 2255, 1897, 1347, 943, 843, 2174, 2460, 1710,
      1015, 841, 815, 1547, 1199, 1238, 951, 2634, 2485, 6871, 6228, 301, 202, 953, 1338, 2115, 2723, 608, 508, 341, 363, 1618, 1627)
  )
  RUnit::checkEquals(
    cx.quality[context=="CHH", sum(as.numeric(pos)), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(273396444558, 282919600073, 137698195394, 168535835281, 191200075125, 212311900157, 75701060631, 81235689888, 74154825642,
      60916711712, 130246222688, 94149778386, 84407645105, 111307951172, 128247086773, 134528051278, 135465135032, 84667838085,
      64917885615, 56758615422, 132389544609, 153880990605, 91566315802, 69198661501, 21001037734, 25962453569, 79021042851,
      72959330970, 71850311932, 56479646308, 84799620854, 94237291846, 159545210605, 164326916078, 12210126684, 6921505496,
      26938878530, 50640456264, 36939107790, 36863076544, 15667795056, 15155859140, 10500573996, 8949271281, 120780623141, 115740654204)
  )
  
  
  # more extended consistency checks on filtered but not thresholded data
  RUnit::checkEquals(
    cx.filt[context=="CG", sum(meth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(140, 151, 67, 61, 107, 116, 83, 55, 69, 0, 44, 21, 136, 207, 74, 40, 86, 90, 4, 29, 63, 49, 140,
      45, 33, 8, 104, 77, 50, 20, 97, 67, 250, 302, 5, 3, 31, 55, 214, 350, 58, 37, 5, 15, 65, 51)
  )
  RUnit::checkEquals(
    cx.filt[context=="CG", sum(unmeth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(503, 656, 263, 284, 393, 625, 135, 92, 91, 74, 361, 266, 231, 299, 298, 273, 299, 84, 308, 192, 403,
      563, 299, 166, 124, 154, 215, 241, 248, 183, 499, 427, 916, 937, 93, 80, 164, 256, 440, 376, 53, 29, 81, 153, 145, 167)
  )
  RUnit::checkEquals(
    cx.filt[context=="CG", sum(as.numeric(pos)), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(70910846257, 61091912169, 28161626105, 25512209140, 36267448482, 46521831699, 28692429340, 19626655012, 22297319654, 8295564267,
      29223919719, 18062482933, 19380130230, 30379246046, 19339970251, 16873851721, 29802187646, 11866520402, 18165011887, 12013972416,
      26086947123, 34450887238, 27850464419, 17035597898, 5856776858, 2574665389, 16902019483, 17605027170, 18984320658, 11683552934,
      22008798733, 20288887142, 37077032663, 29766667492, 3011229920, 3211483594, 5906395308, 11194458997, 11197764917, 10806596119,
      2792182961, 2106023495, 2473611271, 4297344545, 13964534749, 14244322610)
  )
  RUnit::checkEquals(
    cx.filt[context=="CHG", sum(meth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(2, 1, 1, 0, 4, 5, 0, 1, 0, 0, 4, 0, 0, 2, 1, 1, 1, 1, 3, 1, 1, 1, 1, 0, 1, 0, 2, 0, 1, 0, 0, 3, 6, 0, 0, 0, 0, 1, 6, 3, 0, 1, 0, 1, 0, 1)
  )
  RUnit::checkEquals(
    cx.filt[context=="CHG", sum(unmeth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(632, 752, 394, 419, 463, 677, 200, 122, 157, 75, 385, 284, 380, 556, 300, 335, 426, 177, 235, 207, 527, 625, 544, 276, 117, 132, 374,
      368, 340, 176, 680, 628, 1387, 1481, 69, 63, 265, 336, 672, 846, 96, 70, 77, 145, 230, 247)
  )
  RUnit::checkEquals(
    cx.filt[context=="CHG", sum(as.numeric(pos)), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(69949902506, 54891405966, 32707418152, 30646498854, 34053849113, 44872005274, 24685544371, 15825836232, 20180136091, 7761966562, 25932911610,
      18255533534, 20336363297, 32972660373, 14693741851, 18918950210, 35180201144, 11551212580, 13446290110, 12068972433, 29207559112, 33717380334,
      28946953350, 22033587387, 4419157601, 2374139440, 20252265694, 20984684663, 21233835190, 10149252847, 26819841105, 26334824616, 44260232826,
      39929531393, 2118292438, 2436577039, 7506095467, 14216273179, 11827330877, 11496184151, 2672655440, 2242191272, 2104025351, 3710588050,
      16830930354, 15197818899)
  )
  RUnit::checkEquals(
    cx.filt[context=="CHH", sum(meth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(7, 7, 2, 3, 2, 2, 1, 2, 1, 0, 2, 2, 2, 8, 2, 2, 3, 2, 2, 1, 1, 4, 4, 1, 1, 2, 5, 0, 3, 2, 6, 5, 13, 10, 0, 0, 0, 4, 7, 4, 0, 0, 1, 0, 0, 1)
  )
  RUnit::checkEquals(
    cx.filt[context=="CHH", sum(unmeth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(1336, 1559, 750, 784, 873, 1543, 332, 271, 331, 191, 764, 491, 780, 1033, 608, 711, 858, 394, 522, 426, 1177, 1280, 1083, 521, 164, 354, 879,
      733, 640, 348, 1301, 1308, 3039, 3264, 150, 180, 559, 837, 1543, 1527, 121, 141, 167, 248, 422, 565)
  )
  RUnit::checkEquals(
    cx.filt[context=="CHH", sum(as.numeric(pos)), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(137097580990, 127882427697, 59467930734, 59370773028, 64721379521, 92038892552, 43966592857, 33797886822, 45987189090, 19701776306, 46681522970,
      35097035441, 40681798463, 58790306511, 31119462330, 39417630830, 70071852061, 22073312155, 26847386609, 24735110323, 65719444242, 70285659788,
      60082933959, 43298440748, 6095494032, 6115898381, 42718856955, 43640375454, 39988988624, 19748452206, 49331215184, 56432560745, 91693962261,
      93344736990, 5117508821, 6416992151, 14900771359, 35818094672, 29708724094, 23909716952, 3840997396, 4509642318, 5065797755, 6157336608,
      30066085704, 32038175066)
  )
  
  
  
  ### single-end
  
  cx.single <- generateCytosineReport(
    system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"),
    threshold.reads=FALSE, filter.reads=FALSE, report.context="CX", verbose=TRUE
  )
  
  RUnit::checkEquals(
    dim(cx.single),
    c(3236, 6)
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.single$context)[c("CHH", "CHG", "CG")]),
    c(2165, 802, 269)
  )
  
  RUnit::checkEquals(
    c(sum(cx.single$meth), sum(cx.single$unmeth)),
    c(355, 3599)
  )
  
  cx.trim <- generateCytosineReport(
    system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"),
    threshold.reads=FALSE, trim=1, report.context="CX", verbose=FALSE
  )
  cx.notrim <- generateCytosineReport(
    system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"),
    threshold.reads=FALSE, trim=0, report.context="CX", verbose=FALSE
  )
  
  RUnit::checkTrue(
    ! identical(cx.trim, cx.notrim)
  )
  
  RUnit::checkTrue(
    identical(
      data.table::merge.data.table(
        cx.trim[, .(rname, strand, pos, context)],
        cx.trim[, .(rname, strand, pos, context)]
      ),
      data.table::merge.data.table(
        cx.trim[,   .(rname, strand, pos, context)],
        cx.notrim[, .(rname, strand, pos, context)]
      )
    )
  )
  
  
  ### real long-read
  cg.quality  <- generateCytosineReport(
    system.file("extdata", "longread.bam", package="epialleleR"),
    threshold.reads=FALSE, filter.reads=FALSE,
    min.mapq=30, min.baseq=20, min.prob=178, verbose=FALSE
  )
  # extended consistency checks
  RUnit::checkEquals(
    cg.quality[context=="CG", sum(meth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(1303, 1405)
  )
  RUnit::checkEquals(
    cg.quality[context=="CG", sum(unmeth), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(744, 504)
  )
  RUnit::checkEquals(
    cg.quality[context=="CG", sum(as.numeric(pos)), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(20700942290, 18673334762)
  )
  
  
  ### simulated long-read
  output.bam <- tempfile(pattern="output-", fileext=".bam")
  
  # C+m and chebi for other mod
  # modified from htslib/test/base_mods/MM-chebi.sam
  simulateBam(
    flag=c(0),
    seq=c("AGCTCTCCAGAGTCGNACGCCATYCGCGCGCCACCA"),
    pos=1,
    Mm=c("C+m,2,2,1,4,1;C+76792,6,7;N+n,15;"),
    Ml=list(as.integer(c(102,128,153,179,161,187,212,169))),
    output.bam.file=output.bam
  )
  cx.report <- generateCytosineReport(output.bam, filter.reads=FALSE, threshold.reads=FALSE, report.context="CX")
  RUnit::checkEquals(
    cx.report[, .(strand, pos, context, meth, unmeth)],
    data.table::data.table(
      strand=factor("+", levels=c("+","-")),
      pos=as.integer(c(3, 5, 7, 8, 14, 18, 20, 21, 25, 27, 29, 31, 32, 34, 35)),
      context=factor(c(2, 2, 2, 6, 7, 7, 2, 2, 7, 7, 7, 2, 2, 2, 2), levels=1:7,
                     labels=c("NA1", "CHH", "NA3", "NA4", "NA5", "CHG", "CG")),
      meth=  as.integer(c(0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0)),
      unmeth=as.integer(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1))
    )
  )
  cx.report <- generateCytosineReport(output.bam, filter.reads=FALSE, threshold.reads=FALSE, report.context="CX",
                                      min.prob=160, highest.prob=FALSE)
  RUnit::checkEquals(
    cx.report[, .(strand, pos, context, meth, unmeth)],
    data.table::data.table(
      strand=factor("+", levels=c("+","-")),
      pos=as.integer(c(3, 5, 7, 8, 14, 18, 20, 21, 25, 27, 29, 31, 32, 34, 35)),
      context=factor(c(2, 2, 2, 6, 7, 7, 2, 2, 7, 7, 7, 2, 2, 2, 2), levels=1:7,
                     labels=c("NA1", "CHH", "NA3", "NA4", "NA5", "CHG", "CG")),
      meth=  as.integer(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1)),
      unmeth=as.integer(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0))
    )
  )
  
  # C+m and chebi for C+m
  # modified from htslib/test/base_mods/MM-chebi.sam
  simulateBam(
    flag=c(0),
    seq=c("AGCTCTCCAGAGTCGNACGCCATYCGCGCGCCACCA"),
    pos=1,
    Mm=c("C+m,2,2,1,4,1;C+27551,6,7;N+n,15;"),
    Ml=list(as.integer(c(102,128,153,179,161,187,212,169))),
    output.bam.file=output.bam
  )
  cx.report <- generateCytosineReport(output.bam, filter.reads=FALSE, threshold.reads=FALSE, report.context="CX")
  RUnit::checkEquals(
    cx.report[, .(strand, pos, context, meth, unmeth)],
    data.table::data.table(
      strand=factor("+", levels=c("+","-")),
      pos=as.integer(c(3, 5, 7, 8, 14, 18, 20, 21, 25, 27, 29, 31, 32, 34, 35)),
      context=factor(c(2, 2, 2, 6, 7, 7, 2, 2, 7, 7, 7, 2, 2, 2, 2), levels=1:7,
                     labels=c("NA1", "CHH", "NA3", "NA4", "NA5", "CHG", "CG")),
      meth=  as.integer(c(0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1)),
      unmeth=as.integer(c(1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0))
    )
  )
  cx.report <- generateCytosineReport(output.bam, filter.reads=FALSE, threshold.reads=FALSE, report.context="CX",
                                      min.prob=160)
  RUnit::checkEquals(
    cx.report[, .(strand, pos, context, meth, unmeth)],
    data.table::data.table(
      strand=factor("+", levels=c("+","-")),
      pos=as.integer(c(3, 5, 7, 8, 14, 18, 20, 21, 25, 27, 29, 31, 32, 34, 35)),
      context=factor(c(2, 2, 2, 6, 7, 7, 2, 2, 7, 7, 7, 2, 2, 2, 2), levels=1:7,
                     labels=c("NA1", "CHH", "NA3", "NA4", "NA5", "CHG", "CG")),
      meth=  as.integer(c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1)),
      unmeth=as.integer(c(1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0))
    )
  )
  
  # modification on both strands
  # modified from htslib/test/base_mods/MM-double.sam
  simulateBam(
    flag=c(0),
    seq=c("AGGATCTCTAGCGGATCGGCGGGGGATATGCCATAT"),
    pos=1,
    Mm=c("C+m,1,3,0;G-m,0,2,0,4;G+o,4;"),
    Ml=list(as.integer(c(128,153,179,115,141,166,192,102))),
    output.bam.file=output.bam
  )
  cx.report <- generateCytosineReport(output.bam, filter.reads=FALSE, threshold.reads=FALSE, report.context="CX")
  RUnit::checkEquals(
    cx.report[meth>0, .(strand, pos, context)],
    data.table::data.table(
      strand=factor(c("-", "+", "-", "-", "-", "+", "+"), levels=c("+","-")),
      pos=as.integer(c(2, 8, 13, 14, 23, 31, 32)),
      context=factor(c(2, 2, 7, 6, 2, 2, 2), levels=1:7,
                     labels=c("NA1", "CHH", "NA3", "NA4", "NA5", "CHG", "CG"))
    )
  )
  
  #modification on + and -, and on the sequenced and opposite strands
  # modified from htslib/test/base_mods/MM-pileup.sam
  simulateBam(
    flag=c(0, 16),
    seq=c("AGCTCTCCAGAGTCGNACGCCATYCGCGCGCCACCA"),
    pos=1,
    Mm=c("C+m,2,2,1,4,1;C+h,6,7;N+n,15,2;",
         "G-m,0,1,4,1,2;G-h,0,7;N-n,17,2;"),
    Ml=list(as.integer(c(128,153,179,204,230,159,6,215,240)),
            as.integer(c(230,204,179,153,128,6,159,240,215))),
    output.bam.file=output.bam
  )
  cx.report <- generateCytosineReport(output.bam, filter.reads=FALSE, threshold.reads=FALSE, report.context="CX")
  RUnit::checkEquals(
    unname(unlist(cx.report[strand=="-", .(sum(meth), sum(unmeth))])),
    c(0,8)
  )
  RUnit::checkEquals(
    cx.report[strand=="+" & meth>=1, .(pos, context)],
    data.table::data.table(
      pos=as.integer(c(7, 18, 21, 32, 35)),
      context=factor(c(2, 7, 2, 2, 2), levels=1:7,
                     labels=c("NA1", "CHH", "NA3", "NA4", "NA5", "CHG", "CG"))
    )
  )
  
  # modification on '+' and on the sequenced and opposite strands
  # modified from htslib/test/base_mods/MM-orient.sam
  simulateBam(
    flag=c(0),
    seq=c("AGGATCTCTAGCGGATCGGCGGGGGATATGCCATAT"),
    pos=1,
    Mm=c("C+m,2,0,0;G-m,3,1,1;"),
    Ml=list(as.integer(c(128,153,179,128,153,179))),
    output.bam.file=output.bam
  )
  cx.report <- generateCytosineReport(output.bam, filter.reads=FALSE, threshold.reads=FALSE, report.context="CX")
  RUnit::checkEquals(
    dim(cx.report),
    c(20,6)
  )
  RUnit::checkEquals(
    cx.report[context=="CG", .(strand, pos, meth, unmeth)],
    data.table::data.table(
      strand=factor(c("+", "-", "+", "-", "+", "-"), levels=c("+","-")),
      pos=as.integer(c(12, 13, 17, 18, 20, 21)),
      meth=  as.integer(rep.int(1, 6)),
      unmeth=as.integer(rep.int(0, 6))
    )
  )
  
  # modification on '-' and on the sequenced and opposite strands
  # modified from htslib/test/base_mods/MM-orient.sam
  simulateBam(
    flag=16,
    seq=c("AGGATCTCTAGCGGATCGGCGGGGGATATGCCATAT"),
    pos=1,
    Mm=c("C+m,5,1,1;G-m,2,0,0;"),
    Ml=list(as.integer(c(128,153,179,128,153,179))),
    output.bam.file=output.bam
  )
  cx.report <- generateCytosineReport(output.bam, filter.reads=FALSE, threshold.reads=FALSE, report.context="CX")
  RUnit::checkEquals(
    dim(cx.report),
    c(20,6)
  )
  RUnit::checkEquals(
    cx.report[context=="CG", .(strand, pos, meth, unmeth)],
    data.table::data.table(
      strand=factor(c("+", "-", "+", "-", "+", "-"), levels=c("+","-")),
      pos=as.integer(c(12, 13, 17, 18, 20, 21)),
      meth=  as.integer(rep.int(1, 6)),
      unmeth=as.integer(rep.int(0, 6))
    )
  )
  
  # simulateBam(
  #   flag=c(0, 16),
  #   seq=c("AGGATCTCTAGCGGATCGGCGGGGGATATGCCATAT",
  #         "ATATGGCATATCCCCCGCCGATCCGCTAGAGATCCT"),
  #   pos=1,
  #   Mm=c("C+m,1,3,0;", "C+m,1,3,0;",
  #        "G-m,0,0,4,3;", "G-m,0,0,4,3;"),
  #   Ml=list(as.integer(c(128,153,179)), as.integer(c(128,153,179)),
  #           as.integer(c(115,141,166,192)), as.integer(c(115,141,166,192))),
  #   output.bam.file=output.bam
  # )
}
