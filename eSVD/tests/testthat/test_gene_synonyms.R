context("Test gene synonyms")

## convert_synonyms is correct

test_that("convert_synonyms works", {
  vec <- c("asdf", "test")
  res <- convert_synonyms(vec)

  expect_true(length(res) == length(vec))
  expect_true(sum(is.na(res)) == 0)
  expect_true(sum(is.null(res)) == 0)
  expect_true(is.character(res))
})

test_that("convert_synonyms replaces valid genes in-place", {
  vec <- c("Abl2", "ENFL1", "E3B1", "AL023024")
  answer <- c("Abl2", "Chrna4", "Abi1", "Actg1")
  res <- convert_synonyms(vec)

  expect_true(all(res == answer))
})

test_that("convert_synonyms leaves invalid genes in-place", {
  vec <- c("Abl2", "ENFL1", "E3B1", "asdfasdf", "AL023024")
  answer <- c("Abl2", "Chrna4", "Abi1", "asdfasdf", "Actg1")
  res <- convert_synonyms(vec)

  expect_true(all(res == answer))
})
