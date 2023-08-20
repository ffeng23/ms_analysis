#R code functions for running modification analysis


#now let's find for each unique stripped sequence what are the possible modifications 
# with mods[2]
runModStats<-function(dt, output.dir, out.file, m,  quantity.field="Ms1.Area")
{
	m_char<-substr(m,1,1)
	sink(file=here(output.dir,paste0(out.file,"_run2.csv")),append=F)
	cat("Modification stat for ",m ,"\n")
	sink()
	#ss<-unique(df2$Stripped.Sequence)[1]
	#ss<-"FKVATPYSLYVCPEGQNVTLTCR" 
	#ss<-"FKVATPYSLYVCPEGQNVTLTCRLLGPVDK"
	#ss<-"GEVQTCSERRPIRQLTFQDLHLHHGGHQAAQTSHDLAQR"
	results<-data.frame()
	for(ss in unique(dt$Stripped.Sequence))
	{
#		cat("doing peptide sequence:",ss,".....\n")
#		sink(file=here(output.dir,out.file),append=T)
#		cat("'=================\n")
#		cat("Peptide sequence ",ss ,"\n")
#		sink()
		temp<-which(dt$Stripped.Sequence==ss)
		temp_df<-dt[temp,]
		#we need to figure out the unique different modification on different
		# amino acids
		m_indices<-gregexpr(text=temp_df$Modified.Sequence2, pattern=m_char,fixed=T)
		
		#now we need to figure out the stripped sequence index of modified aa
		m_indices_raw<-lapply(m_indices, FUN=
			function(x){
				y=x
				if(x[1]!=-1){
					y=(c(1:length(x))-1)*nchar(m)  # count for each position what is the raw index by subtracting the modified string
					y=x-y-1 # minus 1 because we are count at this point the starting position of "("
				}
				return (y)
			}
		 )

		m_indices_raw_unique<-unique(unlist(m_indices_raw))
		m_indices_raw_unique<-m_indices_raw_unique[-which(m_indices_raw_unique==-1)]	
		m_indices_raw_unique<-sort(m_indices_raw_unique)

		#for each we generate a stats 
		#m_raw=21
		for(m_raw in m_indices_raw_unique)
		{
			
			#get the string with modification *only*
			m_string<-paste0(substr(ss,1,m_raw), m, substr(ss, m_raw+1, nchar(ss)))
			m_string_start<-paste0(substr(ss,1,m_raw), m)
			cat("\t**doing modified sequence:", m_string_start,"\n")
			#get stats
			temp_df2<-temp_df %>% 
				mutate(Modified=grepl(pattern=m_string_start,x=Modified.Sequence2, fixed=T))%>%
				mutate_at(c("Run2","Run","Repeat"), as.factor) 
				
			stats <- temp_df2 %>% 
				#mutate(Modified=grepl(pattern=m_string_start,x=Modified.Sequence2, fixed=T)) %>%
				mutate_at(c("Run2","Run","Repeat", "Modified"), as.factor) %>% 
				group_by(Run2, Repeat, Modified, .drop=F) %>% 
				summarize(Mod_Level=sum(get(quantity.field))
					)
			stats2 <- stats %>% 
				#mutate(Modified=grepl(pattern=m_string_start,x=Modified.Sequence2, fixed=T)) %>%
				#mutate_at(c("Run2","Run","Repeat"), as.factor) %>% 
				group_by(Run2, Repeat) %>% 
				mutate(pct=Mod_Level/sum(Mod_Level)*100, total=sum(Mod_Level))  %>% ungroup %>%
				replace(is.na(.),0)

			stats3 <- stats2 %>% dplyr::filter(Modified==TRUE)%>%
				group_by (Run2) %>% 
				summarize(mean_pct=mean(pct), std_pct=sd(pct)

					) 
			stats3_b<- stats2 %>% dplyr::filter(Modified==TRUE)#%>%
				#group_by (Run2) %>% 
				#summarize(mean_pct=mean(pct), std_pct=sd(pct)


			stats4<-stats3
			stats4$Sequence<-m_string
			results<-rbind(results,stats4)
#			sink(file=here(output.dir,out.file),append=T)
#				cat("'***************\n")
#				cat("\tModification at position ",m_raw ,"\n")
#				cat("\tModified sequence:", m_string,"\n")
#				write.csv(stats2)
#				cat("\tnote: Mod_Level is the sum of Ms1.Area; pct is % of the total.\n")
#				write.csv(stats3)
#				cat("\tnote: mean is the mean % of 3 replicates; and sd is the standard deviation.\n")

#			sink()
			#this last step is very important: remove the modifications 
			# this will make sure the patter in the following step 
			# can be matched.
			temp_df$Modified.Sequence2 <-sub(x=temp_df$Modified.Sequence2, 
					pattern=m_string_start,replacement=substr(ss,1,m_raw), fixed=T)
		} #end of inner loop

	}# end of outer loop
	r.mean<-results %>% pivot_wider(id_cols=Sequence,
		names_from=Run2, values_from=mean,values_fill=0) 
	r.std<-results %>% pivot_wider(id_cols=Sequence,
		names_from=Run2, values_from=std,values_fill=0) 
	
	write.csv(r.mean, file=
			here(output.dir,paste0(out.file,"_mean_reformat.csv"))
		)
	write.csv(r.std, file=
			here(output.dir,paste0(out.file,"_std_reformat.csv"))
			)

}#end of function.


#now let's find for each unique stripped sequence what are the possible modifications 
# with mods[2]
runModStats_counts<-function(dt, output.dir, out.file, m,  quantity.field="Ms1.Area")
{
	m_char<-substr(m,1,1)
	sink(file=here(output.dir,paste0(out.file,"_run2.csv")),append=F)
	cat("Modification stat for ",m ,"\n")
	sink()
	#ss<-unique(df2$Stripped.Sequence)[1]
	#ss<-"FKVATPYSLYVCPEGQNVTLTCR" 
	#ss<-"FKVATPYSLYVCPEGQNVTLTCRLLGPVDK"
	#ss<-"GEVQTCSERRPIRQLTFQDLHLHHGGHQAAQTSHDLAQR"
	results<-data.frame()
	for(ss in unique(dt$Stripped.Sequence))
	{
#		cat("doing peptide sequence:",ss,".....\n")
#		sink(file=here(output.dir,out.file),append=T)
#		cat("'=================\n")
#		cat("Peptide sequence ",ss ,"\n")
#		sink()
		temp<-which(dt$Stripped.Sequence==ss)
		temp_df<-dt[temp,]
		#we need to figure out the unique different modification on different
		# amino acids
		m_indices<-gregexpr(text=temp_df$Modified.Sequence2, pattern=m_char,fixed=T)
		
		#now we need to figure out the stripped sequence index of modified aa
		m_indices_raw<-lapply(m_indices, FUN=
			function(x){
				y=x
				if(x[1]!=-1){
					y=(c(1:length(x))-1)*nchar(m)  # count for each position what is the raw index by subtracting the modified string
					y=x-y-1 # minus 1 because we are count at this point the starting position of "("
				}
				return (y)
			}
		 )

		m_indices_raw_unique<-unique(unlist(m_indices_raw))
		m_indices_raw_unique<-m_indices_raw_unique[-which(m_indices_raw_unique==-1)]	
		m_indices_raw_unique<-sort(m_indices_raw_unique)

		#for each we generate a stats 
		#m_raw=21
		for(m_raw in m_indices_raw_unique)
		{
			
			#get the string with modification *only*
			m_string<-paste0(substr(ss,1,m_raw), m, substr(ss, m_raw+1, nchar(ss)))
			m_string_start<-paste0(substr(ss,1,m_raw), m)
			cat("\t**doing modified sequence:", m_string_start,"\n")
			#get stats
			temp_df2<-temp_df %>% 
				mutate(Modified=grepl(pattern=m_string_start,x=Modified.Sequence2, fixed=T))%>%
				mutate_at(c("Run2","Run","Repeat"), as.factor) 
				
			stats <- temp_df2 %>% 
				#mutate(Modified=grepl(pattern=m_string_start,x=Modified.Sequence2, fixed=T)) %>%
				mutate_at(c("Run2","Run","Repeat", "Modified"), as.factor) %>% 
				group_by(Run2, Repeat, Modified, .drop=F) %>% 
				summarize(Mod_Level=sum(get(quantity.field))
					)
			stats2 <- stats %>% 
				#mutate(Modified=grepl(pattern=m_string_start,x=Modified.Sequence2, fixed=T)) %>%
				#mutate_at(c("Run2","Run","Repeat"), as.factor) %>% 
				group_by(Run2, Repeat) %>% 
				mutate(pct=Mod_Level/sum(Mod_Level)*100, total=sum(Mod_Level))  %>% ungroup %>%
				replace(is.na(.),0)

			stats3 <- stats2 %>% dplyr::filter(Modified==TRUE)%>%
				group_by (Run2) %>% 
				summarize(mean_pct=mean(pct), std_pct=sd(pct)

					) 
			stats3_b<- stats2 %>% dplyr::filter(Modified==TRUE)#%>%
				#group_by (Run2) %>% 
				#summarize(mean_pct=mean(pct), std_pct=sd(pct)


			stats4<-stats3_b
			stats4$Sequence<-m_string
			results<-rbind(results,stats4)
#			sink(file=here(output.dir,out.file),append=T)
#				cat("'***************\n")
#				cat("\tModification at position ",m_raw ,"\n")
#				cat("\tModified sequence:", m_string,"\n")
#				write.csv(stats2)
#				cat("\tnote: Mod_Level is the sum of Ms1.Area; pct is % of the total.\n")
#				write.csv(stats3)
#				cat("\tnote: mean is the mean % of 3 replicates; and sd is the standard deviation.\n")

#			sink()
			#this last step is very important: remove the modifications 
			# this will make sure the patter in the following step 
			# can be matched.
			temp_df$Modified.Sequence2 <-sub(x=temp_df$Modified.Sequence2, 
					pattern=m_string_start,replacement=substr(ss,1,m_raw), fixed=T)
		} #end of inner loop

	}# end of outer loop
	#r.mean<-results %>% pivot_wider(id_cols=Sequence,
	#	names_from=Run2, values_from=mean,values_fill=0) 
	#r.std<-results %>% pivot_wider(id_cols=Sequence,
	#	names_from=Run2, values_from=std,values_fill=0) 
	r.counts<-results %>% pivot_wider(id_cols=Sequence,
		names_from=c(Run2, Repeat), values_from=c(Mod_Level,total),
			values_fill=0)
	write.csv(r.mean, file=
			here(output.dir,paste0(out.file,"_counts_reformat.csv"))
		)
	#write.csv(r.std, file=
	#		here(output.dir,paste0(out.file,"_std_reformat.csv"))
	#		)

}#end of function.