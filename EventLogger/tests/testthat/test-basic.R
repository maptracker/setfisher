library("EventLogger")

test_that("Basic Logging", {
    el <- EventLogger()
    expect_identical(nrow(el$log), 1L, "Empty log setup")
    expect_identical(el$log$Message, "Log initialized")
    expect_message(x <- el$message("hi"))
    expect_identical(el$log$Message[2], "hi")

    

    expect_identical(el$verbose(), TRUE)
    el$verbose(FALSE)
    expect_identical(el$verbose(), FALSE)

    expect_identical(el$useColor(), TRUE)
    el$useColor(FALSE)
    expect_identical(el$useColor(), FALSE)

    expect_identical(el$tidyTime(1/100), "10.000 ms")
    expect_identical(el$tidyTime(1), "1.000 s")
    expect_identical(el$tidyTime(100), "1.667 min")
    
})
