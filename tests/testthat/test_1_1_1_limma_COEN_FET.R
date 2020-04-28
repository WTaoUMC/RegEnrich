context("Initializing, limma, COEN, FET/GSEA, and ranking")

test_that("limma_COEN_FET/GSEA works",{
  data(Lyme_GSE63085)
  data(TFs)
  
  data = log2(Lyme_GSE63085$FPKM + 1)
  x = apply(data, 1, sd)
  data1 = data[x > 0, ]
  data1 = data1[1:2000, ]

  # Lyme_GSE63085$sampleInfo$week = as.factor(Lyme_GSE63085$sampleInfo$week)
  # Lyme_GSE63085$sampleInfo$patientID = as.factor(sub("(\\d+)-(\\d+)", "\\1_\\2",
  #                                                    Lyme_GSE63085$sampleInfo$patientID))

  #### limma ####
  # constructing a RegenrichSet object
  design = model.matrix(~0 + patientID + week, data = Lyme_GSE63085$sampleInfo)
  object = RegenrichSet(expr = data1,
                        colData = Lyme_GSE63085$sampleInfo,
                        method = "limma", minMeanExpr = 0,
                        design = design,
                        contrast = c(rep(0, ncol(design) - 1), 1),
                        networkConstruction = "COEN", # trace = FALSE, nbTrees = 1000, #fast = TRUE,
                        enrichTest = "FET")

  expect_s4_class(object, "RegenrichSet") # test 1

  tmp = capture.output(object <- regenrich_diffExpr(object))
  # expect_equal(object@resDEA@pFC$p[1], 0.8674833214) # test 2
  expect_equal(S4Vectors::mcols(object)$p[1], 0.8674833214) # test 2

  set.seed(1234)
  tmp = capture.output(object <- regenrich_network(object))
  expect_equal(object@topNetwork@elementset$weight[1], 0.03712355302) # test 3

  tmp = capture.output(object <- regenrich_enrich(object, enrichTest = "FET"))
  expect_equal(log(object@resEnrich@allResult$pvalue[1]), log(2.179867472e-09)) # test 4

  tmp = capture.output(object <- regenrich_rankScore(object))
  expect_equal(object@resScore$score[1], 1.679524525) # test 5

  set.seed(1234)
  tmp = capture.output(object <- regenrich_enrich(object, enrichTest = "GSEA"))
  expect_equal(log(round(object@resEnrich@allResult$pval[1], digits = 4)), 
               log(0.0001)) # test 6

  tmp = capture.output(object <- regenrich_rankScore(object))
  expect_equal(object@resScore$score[1], 1.784030808) # test 7
})
