#' Calculate standardized mass index (SMI)
#'
#' Calculates SMI using raw length and mass data
#'
#' @param mass A vector of mass values
#' @param length A vector of length values
#' @param L0 A constant value used in standardization. If not specified, the mean of `length` is used
#' @param resid (Default = FALSE) If `TRUE`, the residuals of ordinary least squares regression also be returned
#'
#' @return A data frame containing the calculated SMI, and optionally, OLS residuals
#' @details This function calculates SMI following Peig and Green. Additionally, the function can return residuals from OLS regression, as an alternative measure / index of body condition. Provide untransformed length and mass values. They will be log-transformed during the calculation of SMI. It is also not necessary to remove observations with missing data values.
#'
#' Peig, J., and A. J. Green. 2009. New perspectives for estimating body condition from mass/length data: The scaled mass index as an alternative method. Oikos 118:1883–1891.
#'
#' @export
#' @examples
#' ## Not run:
#' data("sal_data")
#'
#' smi_results1 <- smi(mass = sal_data$Mass,
#'                     length = sal_data$SVL)
#'
#' smi_results2 <- smi(mass = sal_data$Mass,
#'                     length = sal_data$SVL,
#'                     L0 = mean(sal_data$SVL, na.rm = TRUE))
#'
#' smi_results3 <- smi(mass = sal_data$Mass,
#'                     length = sal_data$SVL,
#'                     L0 = mean(sal_data$SVL, na.rm = TRUE),
#'                     resid = TRUE)
#'
#' plot(smi_results3$smi ~ smi_results3$resid)
#'
#' ## End (Not run)
#'
#' @usage smi(mass, length, L0 = NULL, resid = FALSE)
#'
#' @author Bill Peterman <Peterman.73@@osu.edu>

smi <- function(mass,
                length,
                L0 = NULL,
                resid = FALSE){
  ## Checks
  if(!is.numeric(mass)){
    stop("Mass must be a vector of numeric values")
  }

  if(!is.numeric(length)){
    stop("Length must be a vector of numeric values")
  }

  if(length(mass) != length(length)){
    stop("Mass and Length data have different numbers of observations")
  }

  ## L0
  if(!is.null(L0)){
    if(!is.numeric(L0) | length(L0) > 1){
      stop("L0 must be a single numeric value")
    }
  } else {
    L0 <- mean(length, na.rm = TRUE)
  }

  ## Fit OLS
  ols <- lm(log(mass) ~ log(length),
            na.action = na.exclude)

  p_cor <- cor.test(y = log(mass),
                    x = log(length),
                    method = "pearson")$estimate[[1]]

  b_ols <- ols$coefficients[[2]]
  b_sma <- b_ols / p_cor

  ## Calculate SMI
smi_calc <- mass * (L0 / length)^b_sma

df_out <- data.frame(length = length,
                     mass = mass,
                     smi = smi_calc)

if(isTRUE(resid)){
  df_out$resid <- resid(ols)
}

return(df_out)
}
