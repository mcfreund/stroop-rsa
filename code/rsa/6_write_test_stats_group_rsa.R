## about ----
## 
## calculates basic test statistics (sum of signed ranks) on RSA model fits.

## setup ----

library(here)
library(magrittr)
library(dplyr)
library(data.table)
library(mikeutils)
source(here("code", "strings.R"))
source(here("code", "read_atlases.R"))


models <- c("tdic", "tdi", "tdclust")
measures <- c("coef", "beta", "partr")


## read ----

for (model.i in models) {
  
  stats.subjs <- fread(
    here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_mmp_pearson_residual_glm-tdic.csv"))
  )
  
  stats.subjs <- stats.subjs[is.analysis.group == TRUE, ]  ## EXCLUDE HELD OUT SUBJECTS!
  stats.subjs %<>% full_join(atlas.key$mmp, by = "roi")
    
    
    ## test ----
    
    ## beta
    
    stats.group.beta <- stats.subjs %>%
      group_by(num.roi, param, y) %>%
      summarize(
        v    = wilcox.test(beta, alternative = "greater")$statistic,
        p    = wilcox.test(beta, alternative = "greater")$p.value,
        m    = mean(beta),  ## must be last in this summarize()
      ) %>%
      group_by(param, y) %>%  ## adjust p-value over all parcels
      mutate(p.fdr = p.adjust(p, method = "fdr")) %>%
      ungroup
    
    stats.group.beta %<>% full_join(atlas.key$mmp, by = "num.roi") %>% as.data.table
    
    ## coef
    
    stats.group.coef <- stats.subjs %>%
      group_by(num.roi, param) %>%
      summarize(
        v    = wilcox.test(coef, alternative = "greater")$statistic,
        p    = wilcox.test(coef, alternative = "greater")$p.value,
        m    = mean(coef),  ## must be last in this summarize()
      ) %>%
      group_by(param) %>%  ## adjust p-value over all parcels
      mutate(p.fdr = p.adjust(p, method = "fdr")) %>%
      ungroup
    
    stats.group.coef %<>% full_join(atlas.key$mmp, by = "num.roi") %>% as.data.table
    
    ## partr
    
    stats.group.partr <- stats.subjs %>%
      group_by(num.roi, param) %>%
      summarize(
        v    = wilcox.test(partr, alternative = "greater")$statistic,
        p    = wilcox.test(partr, alternative = "greater")$p.value,
        m    = tanh(mean(atanh(partr))),  ## must be last in this summarize()
      ) %>%
      group_by(param) %>%  ## adjust p-value over all parcels
      mutate(p.fdr = p.adjust(p, method = "fdr")) %>%
      ungroup
    
    stats.group.partr %<>% full_join(atlas.key$mmp, by = "num.roi") %>% as.data.table
    
    
    ## write ----
    
    fwrite(
      stats.group.partr,
      here("out", "rsa", "stats", paste0("group_pro_bias_acc-only_mmp_pearson_residual_glm-", model.i, "_coef.csv"))
    )
    
    fwrite(
      stats.group.partr,
      here("out", "rsa", "stats", paste0("group_pro_bias_acc-only_mmp_pearson_residual_glm-", model.i, "_beta.csv"))
    )
    
    fwrite(
      stats.group.partr,
      here("out", "rsa", "stats", paste0("group_pro_bias_acc-only_mmp_pearson_residual_glm-", model.i, "_partr.csv"))
    )
  
}
  
  
  
  
  
  
}


