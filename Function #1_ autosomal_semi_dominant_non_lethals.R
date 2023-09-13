# This code was written by Ian Hastings, Liverpool School of Tropical Medicine.

# This is one of the functions used in the paper to be published in G3 i.e.
#Hastings, I.M. (2023). "The impact of Fisher’s Reproductive Compensation on raising equilibrium frequencies of semi-dominant, non-lethal mutations under mutation/selection balance".
# Genes, Genomes, Genetics, in press.

#A bit crude but its transparent. If computational efficiency is a priority you will need to alter the code


Find_equilib_freq<-function(q, max_gens, precision, h, s, CR, u ){ 
#NB 'q' is frequency of the mutant allele
#this function will be interrogated by the R function "uniroot" which should rerun the equilibrium frequency of the mutation
#Note: I assume only a single root exists but must revisit and e.g. use "polyroot" or "uniroot.all" to check for multiple roots
#Note:that "q" needs to be the first argument in this function because this is the parameter Uniroot varies to find equilibrium.
#Note: some calculations (e.g. of Size_brood_after_selection) do not depend on q so, in principle,
#these can be calculated before the function is called and passed to the function a vector.
#This would save some computational time as uniroot would not have to recalculate it it at each invocation.
#I just calculate repeatedly in-function on the basis that simplest-is-safest but may revisit is calculations are slow.

#I include a few error checks e.g. check_3_1 which should sum to 1.0, and can be checked by editing the function to catch them
  
  
  
  
#set up a few arrays needed in this function
adult_freq<-c(1,2,3) 
new_adult_freq<-c(1,2,3)
mating_freq<-c(1,2,3,4,5,6) 
brood_genotype<-array(dim=c(6, 3))
Size_brood_after_selection<-c(1,2,3,4,5,6)
Size_brood_after_RC<-c(1,2,3,4,5,6)
freq_change<-c(1,2,3)

lower_limit=1-precision; upper_limit=1+precision # do it once here rather every generation loop


#I am hoping can set up adult freqs according to HW, then return mutant allele frequency
#and that uniroot may still work even if genotypes at end of cycle will not be in HW
#so best to run for e.g. 100 gens to ensure H-W reached
#If uniroot does not work will have to cycle to equilib but this may take longer



#first  calculate all the within brood parameters which do not change with adult_freq so do not re-calculate them every generation

#first index is brood number i.e. 1 to 6
#second index is genotype: 1=++, 2=+m, 3=mm

#mating type #1 (++ with ++)
brood_genotype[1,1]=(1-2*u)
brood_genotype[1,2]=2*u*(1-h*s)  
brood_genotype[1,3]= 0 
check_3_1=brood_genotype[1,1]+brood_genotype[1,2]+brood_genotype[1,3] # set s=1 and this should sum to 1. Checking I have gone the mutations correct
Norm=brood_genotype[1,1]+brood_genotype[1,2]+brood_genotype[1,3]
brood_genotype[1,1]=brood_genotype[1,1]/Norm
brood_genotype[1,2]=brood_genotype[1,2]/Norm
brood_genotype[1,3]=brood_genotype[1,3]/Norm
Size_brood_after_selection[1] = Norm



#mating type #2 (++ with +m)
brood_genotype[2,1]=0.5*(1-2*u)
brood_genotype[2,2]=(0.5*(1-u)+0.5*2*u)*(1-h*s)
brood_genotype[2,3]=0.5*u*(1-s)
check_3_2=brood_genotype[2,1]+brood_genotype[2,2]+brood_genotype[2,3] # set s=1 and this should sum to 1. Checking I have gone the mutations correctNorm<=brood_genotype[2,1]+brood_genotype[2,2]+brood_genotype[2,3]
Norm=brood_genotype[2,1]+brood_genotype[2,2]+brood_genotype[2,3]
brood_genotype[2,1]=brood_genotype[2,1]/Norm
brood_genotype[2,2]=brood_genotype[2,2]/Norm
brood_genotype[2,3]=brood_genotype[2,3]/Norm
Size_brood_after_selection[2] = Norm

#mating type #3 (++ with mm)
brood_genotype[3,1]=0
brood_genotype[3,2]=(1-u)*(1-h*s)
brood_genotype[3,3]=u*(1-s)
check_3_3=brood_genotype[3,1]+brood_genotype[3,2]+brood_genotype[3,3] # set s=1 and this should sum to 1. Checking I have gone the mutations correct
Norm=brood_genotype[3,1]+brood_genotype[3,2]+brood_genotype[3,3]
brood_genotype[3,1]=brood_genotype[3,1]/Norm
brood_genotype[3,2]=brood_genotype[3,2]/Norm
brood_genotype[3,3]=brood_genotype[3,3]/Norm
Size_brood_after_selection[3] = Norm


#mating type #4 (+m with +m)
brood_genotype[4,1]=0.25*(1-2*u)
brood_genotype[4,2]=(0.5*(1-u)+0.25*2*u)*(1-h*s)
brood_genotype[4,3]=(0.25+0.5*u)*(1-s)
check_3_4=brood_genotype[4,1]+brood_genotype[4,2]+brood_genotype[4,3] # set s=1 and this should sum to 1. Checking I have gone the mutations correct
Norm=brood_genotype[4,1]+brood_genotype[4,2]+brood_genotype[4,3]
brood_genotype[4,1]=brood_genotype[4,1]/Norm
brood_genotype[4,2]=brood_genotype[4,2]/Norm
brood_genotype[4,3]=brood_genotype[4,3]/Norm
Size_brood_after_selection[4] = Norm


#mating type #5 (+m with mm)
brood_genotype[5,1]=0
brood_genotype[5,2]=0.5*(1-u)*(1-h*s)
brood_genotype[5,3]=(0.5+0.5*u)*(1-s)
check_3_5=brood_genotype[5,1]+brood_genotype[5,2]+brood_genotype[5,3] # set s=1 and this should sum to 1. Checking I have gone the mutations correct
Norm=brood_genotype[5,1]+brood_genotype[5,2]+brood_genotype[5,3]
brood_genotype[5,1]=brood_genotype[5,1]/Norm
brood_genotype[5,2]=brood_genotype[5,2]/Norm
brood_genotype[5,3]=brood_genotype[5,3]/Norm
Size_brood_after_selection[5]= Norm

#mating type #6 (mm with mm)
brood_genotype[6,1]=0
brood_genotype[6,2]=0
brood_genotype[6,3]=(1-s)
check_3_6=brood_genotype[6,1]+brood_genotype[6,2]+brood_genotype[6,3] # set s=1 and this should sum to 1. Checking I have gone the mutations correct
Norm=brood_genotype[6,1]+brood_genotype[6,2]+brood_genotype[6,3]
brood_genotype[6,1]=brood_genotype[6,1]/Norm
brood_genotype[6,2]=brood_genotype[6,2]/Norm
brood_genotype[6,3]=brood_genotype[6,3]/Norm
Size_brood_after_selection[6] = Norm

#Now to update "Size_brood_after_selection" to allow fRC and check brood size does not exceed 1
temp_vec<-Size_brood_after_selection*CR
Size_brood_after_RC<-replace(temp_vec, temp_vec>1, 1.0) #set max size of brood in the vector to be 1

#calculate frequency of each adult phenotype in the mating population from freq input into the function. Note assumes H-W
#index [1]= freq of ++ genotype, index[2]=freq of +m genotype, index[3]=freq of ++ genotype, 
p<-(1-q)
adult_freq[1]=p^2; adult_freq[2]=2*p*q; adult_freq[3]=q^2; 
check_1a<-sum(adult_freq)



for (ii in 1:max_gens){

#Step1:  find the frequency of each mating type this generation
  check_1b<-sum(adult_freq)
mating_freq[1]<-adult_freq[1]*adult_freq[1]
mating_freq[2]<-2*adult_freq[1]*adult_freq[2]
mating_freq[3]<-2*adult_freq[1]*adult_freq[3]
mating_freq[4]<-adult_freq[2]*adult_freq[2]
mating_freq[5]<-2*adult_freq[2]*adult_freq[3]
mating_freq[6]<-adult_freq[3]*adult_freq[3]
check_2<-sum(mating_freq)


#step 2: calculation frequency of parents of next generation
new_adult_freq[1]=0; new_adult_freq[2]=0; new_adult_freq[3]=0 #recall [1] is ++ [2] is +m [3] is mm

for (mating_type in 1:6) {
new_adult_freq[1]=new_adult_freq[1]+mating_freq[mating_type]*Size_brood_after_RC[mating_type]*brood_genotype[mating_type,1]
new_adult_freq[2]=new_adult_freq[2]+mating_freq[mating_type]*Size_brood_after_RC[mating_type]*brood_genotype[mating_type,2]
new_adult_freq[3]=new_adult_freq[3]+mating_freq[mating_type]*Size_brood_after_RC[mating_type]*brood_genotype[mating_type,3]
}
norm_final=new_adult_freq[1]+new_adult_freq[2]+new_adult_freq[3] #need to normalise because "Size_brood_after_RC" does not sum to unity
new_adult_freq[1]=new_adult_freq[1]/norm_final
new_adult_freq[2]=new_adult_freq[2]/norm_final
new_adult_freq[3]=new_adult_freq[3]/norm_final

check_5<-sum(new_adult_freq)


#now see whether equilib has been reached i.e. frequency change this generation is negligible i.e. less than defined precision
freq_change[1]=new_adult_freq[1]/adult_freq[1]
freq_change[2]=new_adult_freq[2]/adult_freq[2]
freq_change[3]=new_adult_freq[3]/adult_freq[3]



if((min(freq_change)<lower_limit) | (max(freq_change)>upper_limit)) { #i.e. equilib not reached
adult_freq[1]=new_adult_freq[1]
adult_freq[2]=new_adult_freq[2]
adult_freq[3]=new_adult_freq[3]
 } else {  
break }


} # end of cycling through max_gens

mutation_freq=0.5*new_adult_freq[2]+new_adult_freq[3]
#if equilib not reached in max_gens need to report this
##TO DO

return(mutation_freq) #use this option is just invoking the function

#return (mutation_freq-q) #use this option if using uniroot


} #end of function "Find_equilib_freq" 




