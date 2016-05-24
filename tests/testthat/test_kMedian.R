test_data1 <- c(8,7,4,10,19,12,6,1)
test_data2 <- c(20,18,2,9,0,2,3,14)
seq = matrix(c(4,1,4,1,4,3,0,1,2,0,4,0,0,4,4,3,0,0,4,4),1,20)
card = matrix(5,1,20)

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
