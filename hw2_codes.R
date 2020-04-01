lorelogram_v2 <- function(data, id.col = 1, time.col = 4, y.col = 3){
  
  packs <- c("dplyr", "data.table", "ggplot2")
  suppressMessages(packs <- sapply(packs, require, character.only = T))
  data <- data[, c(id.col, time.col, y.col), with = F]
  names(data) <- c("id", "time", "y")
  data <- data.table(data)
  
  # Get the max number of occasion times across subjects
  
  n_rep <- with(data, tapply(time, id, length))
  n_rep <- max(n_rep)
  
  # Get all possible ways of choosing two points ignoring order
  # Note the first column will always have the first time point
  my_comb <- combn(1:n_rep, m = 2)
  my_comb <- t(my_comb)
  
  out_id <- unique(data[, id])
  # For each possible time combination we obtain y_early and y_late across individuals
  out <- lapply(1:nrow(my_comb), function(i){
    
    # Get outcomes at each timepoint 
    y_late <- data[time == my_comb[i, 2], .(id, y)]
    y_early <- data[time == my_comb[i, 1], .(id, y)]
    time_diff <- my_comb[i, 2] - my_comb[i, 1]
    # out_data <- data.table(id = out_id, y_late, y_early, time_diff)
    out_data <- inner_join(data.table(id = out_id), y_early, by = "id") %>% 
      inner_join(data.table(y_late, time_diff),  by = "id") %>% data.table
    names(out_data) <- c("id", "y_late", "y_early", "time_diff")
    out_data
    
    
  })
  out <- do.call(rbind, out)
  out <- out[order(id)]
  
  
  #Predict past outcomes from future by utilizing time differences
  outcome_model <- glm(y_late ~ y_early:factor(time_diff), data=out, family=binomial)
  
  #grab the parameter estimates (ignoring intercept)
  LOR_estimates <- data.table(time_diff = sort(unique(out$time_diff)),
                              point_est = summary(outcome_model)$coef[-1,1])
  
  {
    transparent_legend =  theme(
      legend.background = element_rect(fill ="transparent"),
      legend.key = element_rect(fill = "transparent",
                                color = "transparent")
    )
    
    remove_grid <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    no_x_axis_label <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }
  
  
  # y_lims <- range(LOR_estimates$point_est)
  # y_lims <- round(y_lims, digits = digits)
  # my_breaks <- seq(y_lims[1], y_lims[2], by = .5)
  
  p1 <- ggplot(LOR_estimates, aes(x = time_diff, y = point_est))+
    geom_line(size = 2) +
    geom_point(size = 3, col = "grey") +
    xlab("Lag") + ylab("Log odds ratio") +
    scale_x_continuous(breaks = seq(1, max(LOR_estimates$time_diff), by = .5)) +
    # scale_y_continuous(breaks = my_breaks, limits = y_lims) +
    ggtitle("Lorelogram") +
    theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 18),
          axis.text.y=element_text(size = 18),
          legend.position = c(0.2,0.9),
          legend.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  return(p1)
}


eform <- function(model, level = .95, digits = 3){
  or <- exp(model$coefficients)
  std_err <- summary(model)$coef[,2] 
  or_std_err <- (std_err*or) # Delta Method
  q <- qnorm(1 - (1-level)/2)
  lb <- summary(model)$coef[,1] - q*std_err
  ub <- summary(model)$coef[,1]  + q*std_err
  
  z <- model$coefficients/std_err
  pval <- 2*pnorm(abs(z), lower.tail = F)
  
  out <- cbind(round(or, digits = digits), 
               round(or_std_err, digits = digits), 
               round(z, digits = digits),
               format(pval, digits = 2, scientific = T),
               round(exp(lb), digits = digits), 
               round(exp(ub), digits = digits))
  out <- noquote(out)
  rownames(out) <- names(or)
  colnames(out) <- c("Odds Ratio", "Std. Err", "z", "P>|z|", paste0(level*100, "%", " lower"), paste0(level*100, "%", " upper"))
  if(any(class(model) %in% "geeglm")){
    if(model$std.err ==  "san.se" ){
      colnames(out) <- c("Odds Ratio", "Robust Std. Err", "z", "P>|z|", paste0(level*100, "%", " lower"), paste0(level*100, "%", " upper"))
    }
    
  }
  
  return(out)
  
}


marginal_effects <- function(model, level = .95, model_data, iter = 1e3, bt_sz = 1e3){
  require(data.table)
  model_data$pred_probs <- predict(model, type = "response")
  
  my_margin <- model_data[, .(Estimate = mean(pred_probs)), by = .(laggedy, gender)]
  my_margin <- my_margin[order(laggedy, gender)]
  
  
  # expit <- function(xb){
  #   1/(1+exp(-xb))
  # }
  
  # Use bootstrap
  boot_vars <- sapply(1:nrow(my_margin), function(i){
    
    my_dat <- my_margin[i,]
    boot_dat <- model_data[laggedy == my_dat$laggedy & gender == my_dat$gender, ]
    
    boot_samps <- sapply(seq_len(iter), function(r){
      samp_dt <- sample(x = boot_dat$pred_probs, size = bt_sz, replace = T)
      mean(samp_dt)
    })
    var(boot_samps)*iter
    
  })
  q <- qnorm(1 - (1-level)/2)
  my_margin$SE <- sqrt(boot_vars)
  my_margin$z <- my_margin$Estimate/my_margin$SE
  my_margin[, `:=`(`P-val` = 2*pnorm(abs(z), lower.tail = F),
                   Upper = Estimate + q*SE,
                   Lower = Estimate + q*SE)]
  
  return(my_margin)
}

my_lincom <- function(lin_com = "as.factor(round)2 + as.factor(round)2:gender", model = gee_ft, eform = F, level = .95){
  
  out <- lapply(lin_com, function(x_lin){
    # remove spacing
    x_lin <- gsub(pattern = " ", "", x_lin)
    if(grepl("\\+", x_lin) & !grepl("-", x_lin) ){
      
      get_names <- trimws(x_lin)
      get_names <- unlist(strsplit(get_names, "\\+"))
      get_names <- trimws(get_names)
      
      # Obtain symbols
      # pattern_2_rm <- gsub(pattern = "\\+", "|", x_lin)
      # pattern_wanted <- gsub(pattern_2_rm, "", x_lin)
      pattern_wanted <- gsub("[^\\+|-]", "", x_lin)
      
      # Check if it's the same length as the number of terms in the model. 
      # it shouldn't be
      if(nchar(pattern_wanted) != length(get_names)){
        pattern_wanted <- paste0("+", pattern_wanted)
      }
      
      # Get location of the pluses and the minuses
      minus <- NULL
      plus <- gregexpr("\\+", pattern_wanted)
      plus <- unlist(plus)
      
    }else if(!grepl("\\+", x_lin) & grepl("-", x_lin)){
      get_names <- trimws(x_lin)
      get_names <- unlist(strsplit(get_names, "-"))
      get_names <- trimws(get_names)
      
      # Obtain symbols
      # pattern_2_rm <- gsub(pattern = "-", "|", x_lin)
      # pattern_wanted <- gsub(pattern_2_rm, "", x_lin)
      pattern_wanted <- gsub("[^\\+|-]", "", x_lin)
      
      # Check if it's the same length as the number of terms in the model. 
      # it shouldn't be
      if(nchar(pattern_wanted) != length(get_names)){
        pattern_wanted <- paste0("+", pattern_wanted)
      }
      
      # Get location of the pluses and the minuses
      minus <- gregexpr("-",pattern_wanted)
      minus <- unlist(minus)
      plus <- NULL
      
    }else{
      get_names <- trimws(x_lin)
      get_names <- unlist(strsplit(get_names, "\\+"))
      get_names <- trimws(get_names)
      get_names <- unlist(strsplit(get_names, "-"))
      get_names <- trimws(get_names)
      
      # Obtain symbols
      # pattern_2_rm <- gsub(pattern = "\\+|-", "|", x_lin)
      # pattern_wanted <- gsub(pattern_2_rm, "", x_lin)
      pattern_wanted <- gsub("[^\\+|-]", "", x_lin)
      
      # Check if it's the same length as the number of terms in the model. 
      # it shouldn't be
      if(nchar(pattern_wanted) != length(get_names)){
        pattern_wanted <- paste0("+", pattern_wanted)
      }
      
      # Get location of the pluses and the minuses
      minus <- gregexpr("\\-",pattern_wanted)
      minus <- unlist(minus)
      plus <- gregexpr("\\+", pattern_wanted)
      plus <- unlist(plus)
      
    }
    
    
    
    
    
    
    coef_names <- names(model$coefficients)
    
    contrs <- rep(0, times = length(coef_names))
    contrs[coef_names %in% get_names[minus]] <- -1
    contrs[coef_names %in% get_names[plus]] <- 1
    
    num <- crossprod(contrs, model$coefficients)
    denom <- crossprod(contrs, summary(model)$cov.unscaled) %*% contrs
    
    z <- c(num/sqrt(denom))
    pval <- pnorm(abs(z), lower.tail = F)*2
    
    q <- qnorm(1 - (1-level)/2)
    lb <- num - q*sqrt(denom)
    ub <- num + q*sqrt(denom) 
    
    if(eform){
      res <- c("exp(b)" = exp(num), "Std. Err" = round(sqrt(denom)*exp(num), digits = 4), 
               "z" = z, "P>|z|" = pval, exp(lb), exp(ub)) 
      
      names(res)[5:6] <- c(paste0(level*100, "%", " lower"),
                           paste0(level*100, "%", " upper"))
      res
    }else{
      res <- c("b" = num, "Std. Err" = sqrt(denom), 
               "z" = z, "P>|z|" = pval,lb, ub)
      names(res)[5:6] <- c(paste0(level*100, "%", " lower"),
                           paste0(level*100, "%", " upper"))
      res
    }
    
  })
  
  out <- do.call(rbind, out)
  row.names(out) <- lin_com
  return(out)
  
}