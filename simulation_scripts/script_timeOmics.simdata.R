# RAW DATA
c1 <- c(0, 0.5,1,1.1,1.2,1.8,2.5,5,9)
c3 <-  c(-2,4, 8, 6,4.5,4,3.9, 3, 1)
c2 <- -c1
c4 <- -c3

list(c1,c2,c3,c4)


c1.0 <-  c1
c1.1 <-  c1*1.5
c1.2 <- (c1-0.3)*0.3
c1.3 <- (c1 +0.5)*0.8
c1.4 <- (c1-1)*1.1

c2.0 <-  c2
c2.1 <-  c2*1.5
c2.2 <- (c2-0.3)*0.3
c2.3 <- (c2 +0.5)*0.8
c2.4 <- (c2-1)*1.1

c3.0 <-  c3
c3.1 <-  c3*1.5
c3.2 <- (c3-0.3)*0.3
c3.3 <- (c3 +0.5)*0.8
c3.4 <- (c3-1)*1.1

c4.0 <-  c4
c4.1 <-  c4*1.5
c4.2 <- (c4-0.3)*0.3
c4.3 <- (c4 +0.5)*0.8
c4.4 <- (c4-1)*1.4

# noise 
c0 <- c(0,0.1,0.05,0,0,0.1,0,0.05,0.1) +1
sd(c0)/mean(c0)

data <- list(c1.0,c1.1,c1.2,c1.3,c1.4,c2.0,c2.1,c2.2,c2.3,c2.4,c3.0,c3.1,c3.2,c3.3,c3.4,c4.0,c4.1,c4.2,c4.3,c4.4, c0)
names(data) <- c("c1.0", "c1.1", "c1.2", "c1.3", "c1.4",
                 "c2.0", "c2.1", "c2.2", "c2.3", "c2.4",
                 "c3.0", "c3.1", "c3.2", "c3.3", "c3.4",
                 "c4.0", "c4.1", "c4.2", "c4.3", "c4.4",
                 "c0")
raw.data <- as.data.frame(data)
data.gather <- raw.data %>% rownames_to_column("time") %>%
    mutate(time = as.numeric(time)) %>%
    gather(sample, value, -time)

# SIM DATA
sd <- 0.3
N_Ind <- 5
set.seed(123)

tmp <- data.gather
for(ind in 1:N_Ind){
    vect <- vector(length = nrow(tmp), mode = "numeric")
    for(x in 1:length(vect)){
        vect[x] <- rnorm(1, mean = tmp$value[x], sd = sd)
    }
    name.c <- names(tmp)
    tmp <- data.frame(tmp, vect)
    colnames(tmp) <- c(name.c, LETTERS[ind])
}

sim.data <- tmp %>% dplyr::select(-c(value)) %>%
    gather(ind, value, -c(sample, time))%>%
    mutate(ind = c(paste0(ind, "_", time))) %>% dplyr::select(-time) %>%
    spread(ind, value) %>% column_to_rownames("sample") %>% t

# medelled data
time <- rownames(sim.data) %>% str_split("_") %>% map_chr(~.x[2]) %>% as.numeric()
sampleID <- rownames(sim.data)
lmms.out <- lmms::lmmSpline(data = sim.data, time = time, sampleID = sampleID)

modelled.data <-  as.data.frame(t(lmms.out@predSpline))

timeOmics.simdata <- list(rawdata = raw.data, sim = sim.data,
                          modelled = modelled.data[,-c0])

save(timeOmics.simdata, file = "./data/timeOmics.simdata.RData")
