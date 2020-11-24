###################################################
# subject:data clean for MARC
# email:jingzhy@shanghaitech.edu.cn
# author:Jing, Zhengyu
###################################################

library(stringr)
library(magrittr)
library(tidyverse)
library(xlsx)


#mark every sample, add a marker--------------------------------------------------------------
marc_group <- function(data_ori){   #adds a corresponding sequence number to each sample
  for_out <- data_ori
  for_out$marker <- NA
  for_out$order <- NA
  group_num <- NA
  
  for (index in 1:length(for_out$Name)) {   
    ifelse(str_detect(for_out[index,1],"#"), 
           group_num <- str_extract(for_out[index,1],"\\d+"),
           ifelse(!str_detect(for_out[index,1],"(Name)|(Mapping)"),
                  for_out$marker[[index]] <- group_num,
                  for_out$marker[[index]] <- NA_real_)
               )
  }
  #sets the data type of the column
  for_out %<>% filter(!is.na(marker))
  for_out$NumReads <- round(as.numeric(for_out$NumReads))
  for_out$TPM <- round(as.numeric(for_out$TPM),1)
  for_out <- as_tibble(for_out)
  for_out$Name <- factor(for_out$Name,levels = unique(for_out$Name))
  for_out$mouse_strain <- substring(for_out$Name,1,4)
  
  #add the order of gRNAs
  for(line in unique(for_out$mouse_strain)){
    for_out$order[for_out$mouse_strain == line] <- 1:length(unique(for_out$Name[for_out$mouse_strain == line]))
  }
  for_out$order <- for_out$order
  
  for_out <- for_out %>% select(-Length, -EffectiveLength)

  return(for_out)
}

#Calculate the proportion of each gRNA------------------------------------------------------
percent_marc <- function(data_group, remove_grna_index = c(1,3), remove_gene = "Cd47"){
  new_data_f <- data_group[0,]
  new_data_f$percent <- 1
  new_data_f$percent_remove_some_grna <- 1
  
  for(i in unique(data_group$marker)){
    single_data_f <- data_group %>% filter(marker == i)
    
    #Calculate the proportion of all gRNA
    single_data_f$percent <- single_data_f$NumReads / sum(single_data_f$NumReads)
    #calculate after remove two kinds of gRNA
    s <- sum(single_data_f$NumReads[-c(remove_grna_index, 
                                       which(
                                         str_detect(single_data_f$Name, regex(remove_gene, ignore_case = T))
                                       )
    )]) 
    
    single_data_f$percent_remove_some_grna <- single_data_f$NumReads / s
    single_data_f$percent_remove_some_grna[remove_grna_index] <- NA
    single_data_f$percent_remove_some_grna[str_detect(single_data_f$Name, remove_gene)] <- NA
    
    #combine
    new_data_f <- rbind(new_data_f, single_data_f)
  }
  
  new_data_f <- new_data_f %>% select(-TPM)
  #  new_data_f$logFC <- NA
  new_data_f$marker <- factor(new_data_f$marker, levels = unique(new_data_f$marker))
  new_data_f$besides_ratio <- NA_real_
  new_data_f$ratio_LFC <- NA_real_
  new_data_f$median_reads <- NA_real_
  new_data_f$median_LFC <- NA_real_
  
  return(new_data_f)
}

#add some annotations----------------------------------------------------------------
marc_DataClean_for_plot <- function(data_nor_logfc){
  data_nor_logfc$Name <- str_remove(data_nor_logfc$Name,"\\d+_")
  data_nor_logfc$Name <- factor(data_nor_logfc$Name, levels = unique(data_nor_logfc$Name))
  data_nor_logfc$serial_number <- NA
  
  #deal with 60-mer and 100-mer separately
  sub_1329 <- data_nor_logfc %>% filter(mouse_strain == 1329)
  sub_other <- data_nor_logfc %>% filter(mouse_strain != 1329)
  
  #60-mer first
  if(length(sub_1329$Name) != 0){
    sub_1329$type <- case_when(
      str_detect(sub_1329$Name,"pol") | str_detect(sub_1329$Name,"kpnb") ~ "down",
      str_detect(sub_1329$Name,"p53") | str_detect(sub_1329$Name,"pten") ~ "up",
      str_detect(sub_1329$Name,"cd19") | str_detect(sub_1329$Name,"cd45") ~ "control",
      str_detect(sub_1329$Name,"cd47") ~ "first gRNA")
    
    grna_table <- data.frame(name = as.character(unique(sub_1329$Name)))
    sub_1329$gene <- as.data.frame(str_split(grna_table$name,"-",simplify = TRUE))[[1]] %>% 
      as.character() %>% rep(length(unique(sub_1329$marker)))
    
    #add serial number
    for(gene in unique(sub_1329$gene)){
      if(gene == "cd47"){
        sub_1329$serial_number[sub_1329$gene == gene] <- "a"
      }else{
        sub_1329$serial_number[sub_1329$gene == gene] <- rep(letters[1:10],length(unique(sub_1329$marker)))
      }
    }
  
  
  if(length(sub_other$Name) != 0){
      sub_other$type <- case_when(
      str_detect(sub_other$Name,"NC|GFP") ~ "control",
      TRUE ~ "unknown"
    )
    sub_other$gene <- sub_other$Name %>% as.character()
    sub_other$serial_number <- "a"
    sub_other$gene[str_detect(sub_other$Name,"NC")] <- "NC"
    sub_other$gene[str_detect(sub_other$Name,"GFP")] <- "GFP"
    
    for(i in unique(sub_other$mouse_strain)){
      sub_other$serial_number[str_detect(sub_other$Name,"NC") & sub_other$mouse_strain == i] <- 
        letters[1:length(unique(sub_other$Name[str_detect(sub_other$Name,"NC") & 
                                                 sub_other$mouse_strain == i]))]
      sub_other$serial_number[str_detect(sub_other$Name,"GFP") & sub_other$mouse_strain == i] <- 
        letters[1:length(unique(sub_other$Name[str_detect(sub_other$Name,"GFP") & 
                                                 sub_other$mouse_strain == i]))]
    }
    
  }
  output <- rbind(sub_1329,sub_other) %>% arrange(mouse_strain)
  
  return(output)
}

#add sample information----------------------------------------------------------------------------------------------
data_clean_for_analysis <- function(data_for_plot, sample_path = "./sample_info.xlsx"){
  
  new_data_f <- data_for_plot
  
  sample_info <- read.xlsx(sample_path, sheetIndex = 1, encoding = "UTF-8")
  new_data_f$sample <- sample_info$sample[match(paste0("BZ",new_data_f$marker), sample_info$index_num)]
  new_data_f$sample <- str_replace(new_data_f$sample,"\\s+","_")
  new_data_f$organ <- str_remove(new_data_f$sample,"\\w?\\d+_")
  new_data_f$mouse_index <- str_extract(new_data_f$sample,"\\w?\\d+")
  new_data_f$mouse_index <- factor(new_data_f$mouse_index, levels = unique(new_data_f$mouse_index))
  
  return(new_data_f)
}

#use sample without cas as control group--------------------------------------------------------------------------
gm_mean = function(x, na.rm = TRUE, zero_propagate = T){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero_propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  }else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

normalization_with_median <- 
  function(data, control_l = T, zero = F){
    #add a unique index
    times <- data %>% mutate(temp = paste0(data$marker,data$batch)) %>% group_by(temp) %>% count(temp)
    data$unique_num <- rep(1:length(unique(paste0(data$marker,data$batch))), 
                           times = times$n[match(unique(paste0(data$marker,data$batch)), times$temp)])
    
    new_data <- data[0,]
    for_out <- data[0,]
    #for geometric mean
    for(line in unique(data$mouse_strain)){
      assign(paste0("median",line),
             data %>% filter(mouse_strain == line) %>% 
               group_by(Name) %>% summarise(mean = gm_mean(NumReads, zero_propagate = zero)))
    }
    
    #process each sample
    for(i in unique(data$unique_num)){
      single_data <- data %>% filter(unique_num == i)
      for_median_data <- get(paste0("median", single_data$mouse_strain[[1]])) 
      
      #calculate the size factor
      if(control_l){
        #extract the reads for calculation
        for_median <- unlist(single_data %>% filter(type == "control") %>% select(NumReads))
        names(for_median) <- unlist(single_data %>% filter(type == "control") %>% select(Name))
      }else{
        for_median <- unlist(single_data %>% filter(order != 1) %>% select(NumReads))     
        names(for_median) <- unlist(single_data %>% filter(order != 1) %>% select(Name))
      }
      #size factor!
      size_factor <- median(for_median / for_median_data$mean[match(names(for_median),for_median_data$Name)])
      #normalization
      if(size_factor > 0.01){
        single_data$median_reads <-  
          as.numeric(unlist(
            (single_data %>% select(NumReads)) / 
              (size_factor)
          ))
        
        single_data$median_reads <- round(single_data$median_reads)
      }
      
      new_data <- rbind(new_data, single_data)
    }
    
    #calculate fold change
    #extract the sample without cas
    control_data <- filter(new_data, cas == "no")
    #calculate each sample
    for(i in unique(data$unique_num)){
      single_data <- new_data %>% filter(unique_num == i)
      
      #find the most similar sample for calculation
      control_first <- 
        control_data %>% 
        filter(mouse_strain == unique(single_data$mouse_strain), order == 1) %>% 
        mutate(diff = abs(percent - single_data$percent[[1]])) %>% 
        select(Name, percent, unique_num, sample, diff) %>% 
        arrange(diff)
      
      if(length(control_first[[1]]) == 0){
        for_out <- rbind(for_out, single_data)
        next
      }
      
      six_nearest_data <- control_first$unique_num[1:6]
      within_ten <- control_first$unique_num[control_first$diff < 0.1]
      #select the one has more samples
      if(length(six_nearest_data) >= length(within_ten)){
        control_index <- six_nearest_data
      }else{
        control_index <- within_ten
      }
      
      control_read <- 
        new_data %>% 
        filter(unique_num %in% control_index) %>% 
        group_by(Name) %>% 
        summarise(num = mean(median_reads, na.rm = T))
      #fold change
      single_data$median_LFC <- 
        log2(single_data$median_reads + 1e-20) - 
        log2(control_read$num + 1e-20)
      for_out <- rbind(for_out, single_data)
    }
    return(for_out)
  }
#a <- normalization_with_median(all_data)

normalization_with_beside <- 
  function(data){
    ratio_data <- data[0,]
    for_out <-  data[0,]
    #process each sample
    for(i in unique(data$unique_num)){
      single_data <- data %>% filter(unique_num == i)
      
      #gRNA at both ends is processed separately
      single_data$besides_ratio[2] <- single_data$NumReads[2] / 
        median(single_data$NumReads[3:6] + 1e-20)
      
      single_data$besides_ratio[3] <- single_data$NumReads[3] / 
        median(single_data$NumReads[c(2,4,5,6)] + 1e-20)
      
      single_data$besides_ratio[length(single_data$Name) - 1] <- 
        single_data$NumReads[length(single_data$Name) - 1] / 
        median(single_data$NumReads[length(single_data$Name) - 2:5] + 1e-20)
      
      single_data$besides_ratio[length(single_data$Name) - 2] <- 
        single_data$NumReads[length(single_data$Name)] / 
        median(single_data$NumReads[c(length(single_data$Name) - 1,length(single_data$Name) - 3:5)] + 1e-20)
      
      #process the other gRNAs
      for(n in 3:(length(single_data$Name) -2 )){
        single_data$besides_ratio[n] <- single_data$NumReads[n] / 
          median(single_data$NumReads[c(n-2,n-1,n+1,n+2)] + 1e-20)
      }
      
      ratio_data <- rbind(ratio_data, single_data)
    }
    
    #calculate fold change
    #extract the sample without cas
    control_data <- filter(ratio_data,cas == "no")
      for(i in unique(data$unique_num)){
      single_data <- ratio_data %>% filter(unique_num == i)
      #find the most similar sample for calculation
      control_first <- 
        control_data %>% 
        filter(mouse_strain == unique(single_data$mouse_strain), order == 1) %>% 
        mutate(diff = abs(percent - single_data$percent[[1]])) %>% 
        select(Name, percent, unique_num, sample, diff) %>% 
        arrange(diff)
      
      if(length(control_first[[1]]) == 0){
        for_out <- rbind(for_out, single_data)
        next
      }
      
      six_nearest_data <- control_first$unique_num[1:6]
      within_ten <- control_first$unique_num[control_first$diff < 0.1]
      #select the one has more samples
      if(length(six_nearest_data) >= length(within_ten)){
        control_index <- six_nearest_data
      }else{
        control_index <- within_ten
      }
      
      control_read <- 
        ratio_data %>% 
        filter(unique_num %in% control_index) %>% 
        group_by(Name) %>% 
        summarise(ratio = mean(besides_ratio))
      
      #fold change
      single_data$ratio_LFC <- 
        log2(single_data$besides_ratio + 1e-20) - 
        log2(control_read$ratio + 1e-20)
      
      for_out <- rbind(for_out, single_data)
    }
    return(for_out)
  }
#b <- normalization_with_beside(a)

