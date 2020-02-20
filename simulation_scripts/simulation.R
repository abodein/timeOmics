# raw
c1 <- c(0, 0.5,1,1.1,1.2,1.8,2.5,5,9)
l1 <- smooth.spline(x = c1, spar = 0.3)
p1 <- predict(l1, seq(1,length(c1), length.out = 100))

c3 <-  c(-2,4, 8, 6,4.5,4,3.9, 3, 1)
l3 <- smooth.spline(x = c3, spar = 0.3)
p3 <- predict(l3, seq(1,length(c3), length.out = 100))

c2 <- -c1
l2 <- smooth.spline(x = c2, spar = 0.3)
p2 <- predict(l2, seq(1,length(c2), length.out = 100))


c4 <- -c3
l4 <- smooth.spline(x = c4, spar = 0.3)
p4 <- predict(l4, seq(1,length(c4), length.out = 100))

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

data <- list(c1.0,c1.1,c1.2,c1.3,c1.4,c2.0,c2.1,c2.2,c2.3,c2.4,c3.0,c3.1,c3.2,c3.3,c3.4,c4.0,c4.1,c4.2,c4.3,c4.4)
names(data) <- c("c1.0", "c1.1", "c1.2", "c1.3", "c1.4",
                 "c2.0", "c2.1", "c2.2", "c2.3", "c2.4",
                 "c3.0", "c3.1", "c3.2", "c3.3", "c3.4",
                 "c4.0", "c4.1", "c4.2", "c4.3", "c4.4")
data <- as.data.frame(data)
data.gather <- data %>% rownames_to_column("time") %>%
    mutate(time = as.numeric(time)) %>%
    gather(sample, value, -time)

ggplot(data.gather, aes(time, value, col = sample)) + geom_line() +
    ggtitle("reference profiles")

#save(data, file = "./sim_raw_data.RData")

# Sample (observation)

sd = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5, 2, 3)
N_Ind = 5
EXP = 10
set.seed(123)

final = list()
for(exp in 1:EXP){
    final[[exp]] = list()
    res.sim = list()
    lmms.data = list()
for( S in sd){
    tmp <- data.gather
    for(ind in 1:N_Ind){
        vect <- vector(length = nrow(tmp), mode = "numeric")
        for(x in 1:length(vect)){
            vect[x] <- rnorm(1, mean = tmp$value[x], sd = S)
        }
        name.c <- names(tmp)
        tmp <- data.frame(tmp, vect)
        colnames(tmp) <- c(name.c, LETTERS[ind])
    }
    res.sim[[as.character(S)]] <- tmp %>% mutate(SD = S)
    lmms.data[[as.character(S)]] <- tmp %>% dplyr::select(-c(value)) %>%
        gather(ind, value, -c(sample, time))%>%
        mutate(ind = c(paste0(ind, "_", time))) %>% dplyr::select(-time) %>%
        spread(ind, value) %>% column_to_rownames("sample") %>% t
}
    final[[exp]][["res.sim"]] <- res.sim
    final[[exp]][["lmms.data"]] <- lmms.data
}


set.seed(123)
for(exp in 1:EXP){
    res.sim = list()
    lmms.data = list()
for( S in sd){
    tmp <- data.gather
    for(ind in 1:1){
        vect <- vector(length = nrow(tmp), mode = "numeric")
        for(x in 1:length(vect)){
            vect[x] <- rnorm(1, mean = tmp$value[x], sd = S)
        }
        name.c <- names(tmp)
        tmp <- data.frame(tmp, vect)
        colnames(tmp) <- c(name.c, LETTERS[ind])
    }
    res.sim[[as.character(S)]] <- tmp %>% mutate(SD = S)
    lmms.data[[as.character(S)]] <- tmp %>% dplyr::select(-c(value)) %>%
        gather(ind, value, -c(sample, time))%>%
        mutate(ind = c(paste0(ind, "_", time))) %>% dplyr::select(-time) %>%
        spread(ind, value) %>% column_to_rownames("sample") %>% t 
}
    final[[exp]][["res.sim.1"]] <- res.sim
    final[[exp]][["sim.data.1"]] <- lmms.data
}



final[[1]]$lmms.data$`2` %>% as.data.frame() %>% rownames_to_column("sample") %>% gather(molecule, value, -sample) %>%
    mutate(time = sample %>% str_split("_") %>% map_chr(~.x[2]) %>% as.numeric) %>%
    mutate(ind = sample %>% str_split("_") %>% map_chr(~.x[1])) %>%
    ggplot(aes(time, value)) + geom_point(size = 0.1) + facet_wrap(~ molecule) +
    geom_line(data = data.gather %>% rename(molecule = sample), size = 0.1) + 
    theme_bw() +
    ggtitle("first experiment, simulated data with reference profile (sd = 2)") 
    


for(exp in 1:EXP){
    lmms.res <- list()
    lmms.data <- final[[exp]][["lmms.data"]]
    for(i in names(lmms.data)){
        d <- lmms.data[[i]]
        dt <- rownames(d) %>% str_split("_") %>% map_chr(~.x[2]) %>% as.numeric()
        tmp <- lmms::lmmSpline(data = d, time = dt,
                                         sampleID = rownames(d), deri = FALSE,
                                         basis = "p-spline", numCores = 2)
        lmms.res[[i]] <- t(tmp@predSpline)
    }
    final[[exp]][["lmms.res"]] <- lmms.res
}
get_correspondance_cluster <- function(X){
    # replace cluster label based on occurence
    tmp <- X[,2:3] %>% table %>%
        as.data.frame() %>%
        spread(cluster, Freq) %>%
        column_to_rownames("first_cluster")
    corresp <- apply(X = tmp, FUN = function(x) { colnames(tmp)[which.max(x)[1]]}, MARGIN = 1) %>%
        as.data.frame() %>%
        set_names("new") %>%
        rownames_to_column("old")
    return(corresp)
}


final[[1]]$lmms.data$`2` %>% as.data.frame() %>% rownames_to_column("sample") %>% gather(molecule, value, -sample) %>%
    mutate(time = sample %>% str_split("_") %>% map_chr(~.x[2]) %>% as.numeric) %>%
    mutate(ind = sample %>% str_split("_") %>% map_chr(~.x[1])) %>%
    ggplot(aes(time, value)) + geom_point(size = 0.1) + facet_wrap(~ molecule) +
    geom_line(data = data.gather %>% rename(molecule = sample), size = 0.1) + 
    
    geom_line(data = final[[1]][["lmms.res"]]$`2` %>% as.data.frame() %>% 
                  rownames_to_column("time") %>% gather(molecule, value, -time) %>% mutate(time = as.numeric(time)), color = "red",  size = 0.1) +
    theme_bw() +
    ggtitle("first experiment, simulated data with reference profile (sd = 2) + modelled data(red)") 



# Cumpute accuracy

## With LMMS

accuracy.data.frame <- matrix(ncol = 5, nrow = 100) %>% as.data.frame() %>% purrr::set_names(c("rep1", "noise", "nb_NA", "rep2", "value"))
zz = 0
for(exp in 1:EXP){
    pref.res <- list()
    lmms.res <- final[[exp]][["lmms.res"]]
    conf.mat <- list()
    for(i in names(lmms.res)){
        zz = zz + 1
        pca.res <- mixOmics::pca(lmms.res[[i]], ncomp = 2)
        pca.res.cluster <- getCluster(pca.res) %>%
            dplyr::select(molecule, cluster) %>%
            mutate(first_cluster = molecule %>% str_split("\\.") %>% map_chr(~.x[1]))
        correspondace_cluster <- get_correspondance_cluster(pca.res.cluster)
        classify.res <- pca.res.cluster %>% left_join(correspondace_cluster, by = c("first_cluster"="old")) %>%
            dplyr::select(cluster, new) %>% lapply( function(x) x %>% as.character() %>% as.numeric())
        classify.res <- lapply(classify.res, function(x) factor(x, levels = sort(unique(classify.res$cluster))))
    
        tt <- table(classify.res) %>% as.data.frame() %>%
            spread(cluster, Freq) %>% arrange(as.numeric(as.character(new))) %>%
            column_to_rownames("new")
        tt <- tt %>% dplyr::select(colnames(tt) %>% as.numeric() %>% sort %>% as.character()) %>% as.matrix()
        conf.mat[[i]] <- tt
        pref.res[[i]] <- sum(diag(tt) / sum(tt))
    }
    final[[exp]][["pref.res"]] <- pref.res
    final[[exp]][["conf"]] <- conf.mat
}

zz = 0
# "rep1", "noise", "nb_NA", "rep2", "value" ## rep1 = rep2 to compare with LMMS inputation
for(exp in 1:EXP){
    for(noise in names(final[[exp]][["pref.res"]])){
        zz = zz + 1
        accuracy.data.frame[zz, ] <- c(exp, noise, 0, exp, final[[exp]][["pref.res"]][[noise]])
    }
}



### Confusion matrix (LMMS)

mean (10 repetitions) confusion matrix for each SD.



for(sd in  names(final[[1]][["conf"]])){
    mat <- matrix(0, ncol = 4, nrow = 4, 
                  dimnames = list(c(-2,-1, 1, 2),c(-2,-1, 1, 2)))
    for( exp in 1:length(final)){
        mat <- mat + final[[exp]][["conf"]][[sd]]
    }
    mat <- mat/length(final[[exp]][["conf"]])
    print(paste0("SD = ", sd))
    print(mat)
}

## Without LMMS

```{r}
for(exp in 1:EXP){
    pref.res <- list()
    conf.mat <- list()
    sim.data <-  final[[exp]][["sim.data.1"]]
    for(i in names(sim.data)){
        pca.res <- mixOmics::pca(sim.data[[i]], ncomp = 2)
        pca.res.cluster <- getCluster(pca.res) %>%
            dplyr::select(molecule, cluster) %>%
            mutate(first_cluster = molecule %>% str_split("\\.") %>% map_chr(~.x[1]))
        correspondace_cluster <- get_correspondance_cluster(pca.res.cluster)
        classify.res <- pca.res.cluster %>% left_join(correspondace_cluster, by = c("first_cluster"="old")) %>%
            dplyr::select(cluster, new) %>% lapply( function(x) x %>% as.character() %>% as.numeric())
        classify.res <- lapply(classify.res, function(x) factor(x, levels = sort(unique(classify.res$cluster))))
    
        tt <- table(classify.res) %>% as.data.frame() %>%
            spread(cluster, Freq) %>% arrange(as.numeric(as.character(new))) %>%
            column_to_rownames("new")
        tt <- tt %>% dplyr::select(colnames(tt) %>% as.numeric() %>% sort %>% as.character()) %>% as.matrix()
        conf.mat[[i]] <- tt
        pref.res[[i]] <- sum(diag(tt) / sum(tt))
    }
    final[[exp]][["pref.res.1"]] <- pref.res
    final[[exp]][["conf.1"]] <- conf.mat
}



### Confusion matrix (no LMMS)

mean (10 repetitions) confusion matrix for each SD.

```{r}

for(sd in  names(final[[1]][["conf.1"]])){
    mat <- matrix(0, ncol = 4, nrow = 4, 
                  dimnames = list(c(-2,-1, 1, 2),c(-2,-1, 1, 2)))
    for( exp in 1:length(final)){
        mat <- mat + final[[exp]][["conf.1"]][[sd]]
    }
    mat <- mat/length(final[[exp]][["conf.1"]])
    print(paste0("SD = ", sd))
    print(mat)
}
```




```{r plot_accuracy}
result <- lapply(final, function(x) x$pref.res) %>% unlist 
result <- as.data.frame(result) %>% set_names("CCV") %>% 
    mutate("SD"=names(result)) %>% group_by(SD) %>% 
    summarise(mu = mean(CCV), error = sd(CCV))

ggplot(result, aes(x = SD, y = mu)) + 
    geom_bar(stat="identity", color="black",  position=position_dodge(), fill = "white", width = 0.5) +
    geom_errorbar(aes(ymin=mu-error, ymax=mu+error), width=.2,
                 position=position_dodge(.9))  +
    ylab("CCV") + ggtitle("Overall accuracy of clustering by noise") +
    theme_bw()
```

## Example of clustering results

### raw data clustering 

Unscaled profiles => Good

```{r}
pca.res <- mixOmics::pca(data, ncomp = 2) 
plotLong(pca.res, time = 1:9, scale = F, center = F)
```

### clustering with noise

not as good

unscale modelled data clustering with SD = 2.

```{r}
pca.res <- mixOmics::pca(final[[1]]$lmms.res$`2`, ncomp = 2)  
plotLong(pca.res, time = 1:9, scale = F, center = F)
```


## Accuracy no LMMS


result.1 <- lapply(final, function(x) x$pref.res.1) %>% unlist 
result.1 <- as.data.frame(result.1) %>% set_names("CCV") %>% 
    mutate("SD"=names(result.1)) %>% group_by(SD) %>% 
    summarise(mu = mean(CCV), error = sd(CCV))

res <- rbind(result %>% mutate(modelling = "LMMS"), result.1 %>% mutate(modelling = "no LMMS")) 

ggplot(res, aes(x = SD, y = mu, fill = modelling)) + 
    geom_bar(stat="identity", color="black",  position=position_dodge(), width = 0.5) +
    geom_errorbar(aes(ymin=mu-error, ymax=mu+error), width=.2,
                 position=position_dodge(.5))  +
    ylab("CCV") + ggtitle("Overall accuracy of clustering by noise") +
    theme_bw()






