# This code was written by Ian Hastings, Liverpool School of Tropical Medicine.

# This is one of the functions used in the paper to be published in G3 i.e.
#Hastings, I.M. (2023). "The impact of Fisher’s Reproductive Compensation on raising equilibrium frequencies of semi-dominant, non-lethal mutations under mutation/selection balance".
# Genes, Genomes, Genetics, in press.

#A bit crude but its transparent. If computational efficiency is a priority you will need to alter the code


Find_dimorphic_sex_linked_equilib_freq<-function(q, max_gens, precision, h, s, CR, u ){ 
#NB 'q' is frequency of the mutant allele


#set up a few arrays needed in this function
adult_male_freq<-c(1,2) 
adult_female_freq<-c(1,2,3)
new_adult_male_freq<-c(1,2)
new_adult_female_freq<-c(1,2,3)
mating_freq<-c(1,2,3,4,5,6) 
brood_male_genotype<-array(dim=c(6, 2))
brood_female_genotype<-array(dim=c(6, 3))
Size_male_brood_after_selection<-c(1,2,3,4,5,6)
Size_female_brood_after_selection<-c(1,2,3,4,5,6)
#Size_female_brood_after_RC<-c(1,2,3,4,5,6)
freq_change<-c(1,2,3,4,5)
equilib_reached=0
lower_limit=1-precision; upper_limit=1+precision # do it once here rather every generation loop

return_data<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
names(return_data)  <- c("equilib_freq", "generations_run", "error_code", #[indices 1-3]
                         "check_1a", "check_1b", "check_2",
                         "check_3_1", "check_3_2", "check_3_3", "check_3_4", "check_3_5", "check_3_6", 
                         "check_5a", "check_5b")

#the function will scan the internal checks to check they sum to unity (except 3).
#If a check fails then will change return_data$error_code from 0 to 1
#they won't sum to exactly 1 due to rounding errors so need to set allowable limits
return_data$error_code=0 
min_lim=0.9999999; max_lim =1.000001          


#first  calculate all the within brood parameters which do not change with adult_freq so do not do them every generation

#first index is brood number i.e. 1 to 6
#second index is genotype: 1=++, 2=+m, 3=mm

#mating type #1 (+Y father with ++ mother)

brood_male_genotype[1,1]=0.5*(1-u)
brood_male_genotype[1,2]=0.5*u*(1-s)

brood_female_genotype[1,1]= 0.5*(1-2*u) 
brood_female_genotype[1,2]=0.5*2*u*(1-h*s)
brood_female_genotype[1,3]=0

norm_male=brood_male_genotype[1,1]+brood_male_genotype[1,2]
norm_female= brood_female_genotype[1,1]+brood_female_genotype[1,2]+brood_female_genotype[1,3]
check_3_1=norm_male+norm_female

brood_male_genotype[1,1]=brood_male_genotype[1,1]/norm_male
brood_male_genotype[1,2]=brood_male_genotype[1,2]/norm_male

brood_female_genotype[1,1]=brood_female_genotype[1,1]/norm_female
brood_female_genotype[1,2]=brood_female_genotype[1,2]/norm_female
brood_female_genotype[1,3]=brood_female_genotype[1,3]/norm_female

Size_male_brood_after_selection[1] = norm_male
Size_female_brood_after_selection[1] = norm_female

#mating type #2 (+Y father with +m mother)
brood_male_genotype[2,1]=0.25*(1-u)
brood_male_genotype[2,2]=(0.25+0.25*u)*(1-s)
  
brood_female_genotype[2,1]=0.25*(1-2*u)  
brood_female_genotype[2,2]=(0.25*(1-u)+0.25*2*u)*(1-h*s)
brood_female_genotype[2,3]=0.25*u*(1-s)
  
norm_male=brood_male_genotype[2,1]+brood_male_genotype[2,2]
norm_female=brood_female_genotype[2,1]+brood_female_genotype[2,2]+brood_female_genotype[2,3]
check_3_2=norm_male+norm_female


brood_male_genotype[2,1]=brood_male_genotype[2,1]/norm_male
brood_male_genotype[2,2]=brood_male_genotype[2,2]/norm_male

brood_female_genotype[2,1]=brood_female_genotype[2,1]/norm_female
brood_female_genotype[2,2]=brood_female_genotype[2,2]/norm_female
brood_female_genotype[2,3]=brood_female_genotype[2,3]/norm_female

Size_male_brood_after_selection[2] = norm_male
Size_female_brood_after_selection[2] = norm_female

#mating type #3 (+Y father with mm mother)
brood_male_genotype[3,1]=0
brood_male_genotype[3,2]=0.5*(1-s)
  
brood_female_genotype[3,1]= 0 
brood_female_genotype[3,2]=0.5*(1-u)*(1-h*s)
brood_female_genotype[3,3]=0.5*u*(1-s)
  
norm_male=brood_male_genotype[3,1]+brood_male_genotype[3,2]
norm_female=brood_female_genotype[3,1]+brood_female_genotype[3,2]+brood_female_genotype[3,3]
check_3_3=norm_male+norm_female


brood_male_genotype[3,1]=brood_male_genotype[3,1]/norm_male
brood_male_genotype[3,2]=brood_male_genotype[3,2]/norm_male

brood_female_genotype[3,1]=brood_female_genotype[3,1]/norm_female
brood_female_genotype[3,2]=brood_female_genotype[3,2]/norm_female
brood_female_genotype[3,3]=brood_female_genotype[3,3]/norm_female

Size_male_brood_after_selection[3] = norm_male
Size_female_brood_after_selection[3] = norm_female

#mating type #4 (mY father with ++ mother)
brood_male_genotype[4,1]=0.5*(1-u)
brood_male_genotype[4,2]=0.5*u*(1-s)
  
brood_female_genotype[4,1]=0
brood_female_genotype[4,2]=0.5*(1-u)*(1-h*s)
brood_female_genotype[4,3]=0.5*u*(1-s)
  
norm_male=brood_male_genotype[4,1]+brood_male_genotype[4,2]
norm_female=brood_female_genotype[4,1]+brood_female_genotype[4,2]+brood_female_genotype[4,3]
check_3_4=norm_male+norm_female

brood_male_genotype[4,1]=brood_male_genotype[4,1]/norm_male
brood_male_genotype[4,2]=brood_male_genotype[4,2]/norm_male

brood_female_genotype[4,1]=brood_female_genotype[4,1]/norm_female
brood_female_genotype[4,2]=brood_female_genotype[4,2]/norm_female
brood_female_genotype[4,3]=brood_female_genotype[4,3]/norm_female

Size_male_brood_after_selection[4] = norm_male
Size_female_brood_after_selection[4] = norm_female

#mating type #5 (mY father with +m mother)
brood_male_genotype[5,1]=0.25*(1-u)
brood_male_genotype[5,2]=(0.25+0.25*u)*(1-s)
  
brood_female_genotype[5,1]=0  
brood_female_genotype[5,2]=0.25*(1-u)*(1-h*s)
brood_female_genotype[5,3]=(0.25+0.25*u)*(1-s)
  
norm_male=brood_male_genotype[5,1]+brood_male_genotype[5,2]
norm_female=brood_female_genotype[5,1]+brood_female_genotype[5,2]+brood_female_genotype[5,3]
check_3_5=norm_male+norm_female


brood_male_genotype[5,1]=brood_male_genotype[5,1]/norm_male
brood_male_genotype[5,2]=brood_male_genotype[5,2]/norm_male

brood_female_genotype[5,1]=brood_female_genotype[5,1]/norm_female
brood_female_genotype[5,2]=brood_female_genotype[5,2]/norm_female
brood_female_genotype[5,3]=brood_female_genotype[5,3]/norm_female

Size_male_brood_after_selection[5] = norm_male
Size_female_brood_after_selection[5] = norm_female

#mating type #6 (mY father with mm mother)
brood_male_genotype[6,1]=0
brood_male_genotype[6,2]=0.5*(1-s)
  
brood_female_genotype[6,1]= 0 
brood_female_genotype[6,2]=0
brood_female_genotype[6,3]=0.5*(1-s)
  
norm_male=brood_male_genotype[6,1]+brood_male_genotype[6,2]
norm_female=brood_female_genotype[6,1]+brood_female_genotype[6,2]+brood_female_genotype[6,3]
check_3_6=norm_male+norm_female # set s=0 and this should sum to 1. Checking I have gone the mutations correct

brood_male_genotype[6,1]=brood_male_genotype[6,1]/norm_male
brood_male_genotype[6,2]=brood_male_genotype[6,2]/norm_male

brood_female_genotype[6,1]=brood_female_genotype[6,1]/norm_female
brood_female_genotype[6,2]=brood_female_genotype[6,2]/norm_female
brood_female_genotype[6,3]=brood_female_genotype[6,3]/norm_female

Size_male_brood_after_selection[6] = norm_male
Size_female_brood_after_selection[6] = norm_female

#Now to update "Size_brood_after_selection" to allow fRC and check brood size does not exceed 1
temp_male_vec<-Size_male_brood_after_selection*CR
Size_male_brood_after_RC<-replace(temp_male_vec, temp_male_vec>0.5, 0.5) #set max size of male brood in the vector to be 0.5

temp_female_vec<-Size_female_brood_after_selection*CR
Size_female_brood_after_RC<-replace(temp_female_vec, temp_female_vec>0.5, 0.5) #set max size of male brood in the vector to be 0.5



#calculate frequency of each adult genotype in the mating population from freq input into the function. Note assumes H-W
#index [1]= freq of ++ genotype, index[2]=freq of +m genotype, index[3]=freq of ++ genotype, 
p<-(1-q)
adult_male_freq[1]=p;adult_male_freq[2]=q;
adult_female_freq[1]=p^2; adult_female_freq[2]=2*p*q; adult_female_freq[3]=q^2; 
check_1a<-sum(adult_male_freq)
check_1b<-sum(adult_female_freq)


for (ii in 1:max_gens){

#Step1:  find the frequency of each mating type this generation

mating_freq[1]<-adult_male_freq[1]*adult_female_freq[1]
mating_freq[2]<-adult_male_freq[1]*adult_female_freq[2]
mating_freq[3]<-adult_male_freq[1]*adult_female_freq[3]
mating_freq[4]<-adult_male_freq[2]*adult_female_freq[1]
mating_freq[5]<-adult_male_freq[2]*adult_female_freq[2]
mating_freq[6]<-adult_male_freq[2]*adult_female_freq[3]
check_2<-sum(mating_freq)


#step 2: calculation frequency of parents of next generation
new_adult_male_freq[1]=0; new_adult_male_freq[2]=0; 
new_adult_female_freq[1]=0;new_adult_female_freq[2]=0;new_adult_female_freq[3]=0;  #recall [1] is ++ [2] is +m [3] is mm

for (mating_type in 1:6) {
new_adult_male_freq[1]=new_adult_male_freq[1]+mating_freq[mating_type]*Size_male_brood_after_RC[mating_type]*brood_male_genotype[mating_type,1]
new_adult_male_freq[2]=new_adult_male_freq[2]+mating_freq[mating_type]*Size_male_brood_after_RC[mating_type]*brood_male_genotype[mating_type,2]

new_adult_female_freq[1]=new_adult_female_freq[1]+mating_freq[mating_type]*Size_female_brood_after_RC[mating_type]*brood_female_genotype[mating_type,1]
new_adult_female_freq[2]=new_adult_female_freq[2]+mating_freq[mating_type]*Size_female_brood_after_RC[mating_type]*brood_female_genotype[mating_type,2]
new_adult_female_freq[3]=new_adult_female_freq[3]+mating_freq[mating_type]*Size_female_brood_after_RC[mating_type]*brood_female_genotype[mating_type,3]
}

norm_male=new_adult_male_freq[1]+new_adult_male_freq[2] #need to normalise because "Size_brood_after_RC" does not sum to unity
norm_female=new_adult_female_freq[1]+new_adult_female_freq[2]+new_adult_female_freq[3]

new_adult_male_freq[1]=new_adult_male_freq[1]/norm_male
new_adult_male_freq[2]=new_adult_male_freq[2]/norm_male

new_adult_female_freq[1]=new_adult_female_freq[1]/norm_female
new_adult_female_freq[2]=new_adult_female_freq[2]/norm_female
new_adult_female_freq[3]=new_adult_female_freq[3]/norm_female

check_5a<-sum(new_adult_male_freq)
check_5b<-sum(new_adult_female_freq)


freq_change[1]=new_adult_female_freq[1]/adult_female_freq[1]
freq_change[2]=new_adult_female_freq[2]/adult_female_freq[2]
freq_change[3]=new_adult_female_freq[3]/adult_female_freq[3]
freq_change[4]=new_adult_male_freq[1]/adult_male_freq[1]
freq_change[5]=new_adult_male_freq[2]/adult_male_freq[2]

#NB checks 1 and 3 done before cycling generations so no need to check each loop
if((check_2> max_lim) | (check_2<min_lim))
  return_data$error_code <- 1
if((check_5a> max_lim) | (check_5a<min_lim))
  return_data$error_code <- 1
if((check_5b> max_lim) | (check_5b<min_lim))
  return_data$error_code <- 1


if(min(freq_change)>lower_limit && max(freq_change)<upper_limit)
{
  equilib_reached=1 
} 

adult_male_freq[1]=new_adult_male_freq[1]
adult_male_freq[2]=new_adult_male_freq[2]

adult_female_freq[1]=new_adult_female_freq[1]
adult_female_freq[2]=new_adult_female_freq[2]
adult_female_freq[3]=new_adult_female_freq[3]

if(equilib_reached==1){
break  #break out out of the loop cycling generations
}

} # end of cycling through max_gens

# do this at the end  return_data$equilib_freq <- 
return_data$generations_run <- ii
#return_data$error_code already stored at this location 
return_data$check_1a <- check_1a
return_data$check_1b <- check_1b
return_data$check_2 <- check_2
return_data$check_3_1 <- check_3_1  #NB the check_3 should only sum to zero when s=0 (i.e. no selection). It checks the mutaion rate are correct
return_data$check_3_2 <- check_3_2
return_data$check_3_3 <- check_3_3
return_data$check_3_4 <- check_3_4
return_data$check_3_5 <- check_3_5
return_data$check_3_6 <- check_3_6
return_data$check_5a <- check_5a
return_data$check_5b <- check_5b

male_mut_freq_adult=adult_male_freq[2]
female_mut_freq_adult= adult_female_freq[2]*0.5+adult_female_freq[3]
mutation_freq_adults=(male_mut_freq_adult/3)+(2*female_mut_freq_adult/3)
mutation_freq_gametes=mutation_freq_adults+(1-mutation_freq_adults)*u

return_data$equilib_freq <-mutation_freq_gametes 

return (return_data) 


} #end of function "Find_dimorphic_sex_linked_equilib_freq" 


