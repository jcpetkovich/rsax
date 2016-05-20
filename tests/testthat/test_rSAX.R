test_data1 <- c(8,7,4,10,19,12,6,1)
test_data2 <- c(20,18,2,9,0,2,3,14)
seq = matrix(c(4,1,4,1,4,3,0,1,2,0,4,0,0,4,4,3,0,0,4,4),1,20)
card = matrix(5,1,20)

test_that("Alphabet Size 4 Output -- iSAX = false", {
  #Segment Size 2, Alphabet Size 4, iSAX = false
  expect_equal(runSAX(test_data1,2,4,iSAX = F), c(1,1,3,0))
})

test_that("Alphabet Size 8 Output -- iSAX = true",{
  #Segment Size 2, Alphabet Size 8, iSAX = true
  expect_equal(runSAX(test_data1,2,8,iSAX = T), list(SAXWord = c(3,3,7,1),Cardinality = c(8,8,8,8)))
})

test_that("Normalize Data Function",{
  expect_equal(runNormData(test_data1),c(-0.0732143, -0.2684524357, -0.8541668653, 0.3172619641, 2.0744052, 0.7077382207, -0.4636905789, -1.4398812), tolerance = 0.0000001)
})

test_that("Conversion to PAA Function",{
  norm_data <- c(-0.0732143, -0.2684524357, -0.8541668653, 0.3172619641, 2.0744052, 0.7077382207, -0.4636905789, -1.4398812)
  expect_equal(runToPAA(norm_data,2),c(-0.1708333641, -0.2684524655, 1.3910716772, -0.9517859221))
})

test_that("Conversion to SAX Function",{
  PAA_data <- c(-0.1708333641, -0.2684524655, 1.3910716772, -0.9517859221)
  expect_equal(runToSAX(PAA_data,8), c(3,3,7,1))
})

test_that("Euclidean Distance Function",{
  expect_equal(eucDis(test_data1,test_data2), 30.14962769)
})

test_that("Minimum Distance Function -- Same Cardinality",{
  expect_equal(minDis(c(3,3,7,1),c(8,8,8,8),c(7,2,1,3),c(8,8,8,8),8), 3.091915369)
})

test_that("Minimum Distance Function -- Different Cardinality",{
  expect_equal(minDis(c(3,3,7,1),c(8,8,8,8),c(0,1,1,0),c(2,2,2,2),8), 0)
  expect_equal(minDis(c(1,2,3,4),c(7,7,7,7),c(1,2,2,1),c(4,4,4,4),8,TRUE), 0.7200495)
})

test_that("kMedian Function -- No Previous Centroids",{
  result1 = list(Sequence = c(2,0,4,0),Cardinality = c(5,5,5,5))
  result2 = list(Sequence = c(0,0,4,4),Cardinality = c(5,5,5,5))
  result3 = list(Sequence = c(4,3,0,0),Cardinality = c(5,5,5,5))
  result = list(result1,result2,result3)
  expect_equal(runKMedian(seq,card,2,4,3,60, matrix(1,1,1),matrix(1,1,1),F), result)
})

test_that("kMedian Function -- With Previous Centroids",{
  prevCent = matrix(c(3,0,2,1,1,4,2,1,0,2,2,3),3,4)
  result1 = list(Sequence = c(2,0,4,0),Cardinality = c(5,5,5,5))
  result2 = list(Sequence = c(0,0,4,4),Cardinality = c(5,5,5,5))
  result3 = list(Sequence = c(4,3,0,0),Cardinality = c(5,5,5,5))
  result = list(result1,result2,result3)
  expect_equal(runKMedian(seq,card,2,4,3,60, prevCent,matrix(5,3,4),T), result)
})
