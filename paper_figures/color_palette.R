library(ggsci)

col.pal <- pal_jama()(7)
genetics.include <- c("FUS-ALS", "sporadic", "Healthy", "ALS")
names(col.pal)[c(4, 5, 1, 3)] <- genetics.include
