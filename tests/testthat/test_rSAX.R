test_data <- c(8,7,4,10,19,12,6,1)

test_that("Alphabet Size 4 Output -- iSAX = false", {
  #Segment Size 2, Alphabet Size 4, iSAX = false
  expect_equal(runSAX(test_data,2,4,iSAX = F), c(1,1,3,0))
})

test_that("Alphabet Size 8 Output -- iSAX = true",{
  #Segment Size 2, Alphabet Size 8, iSAX = true
  expect_equal(runSAX(test_data,2,8,iSAX = T), list(SAXWord = c(3,3,7,1),Cardinality = c(8,8,8,8)))
})

test_that("Normalize Data Function",{
  expect_equal(runNormData(test_data),c(-0.0732143, -0.2684524357, -0.8541668653, 0.3172619641, 2.0744052, 0.7077382207, -0.4636905789, -1.4398812))
})

test_that("Conversion to PAA Function",{
  expect_equal(runToPAA(test_data,2),c(-0.1708333641, -0.2684524655, 1.3910716772, -0.9517859221))
})

