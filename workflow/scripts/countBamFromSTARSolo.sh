MAQ=${MAQ:="255"}
NUMHIT=${NUMHIT:=10000}
SORTCPU=${SORTCPU:=8}
mawk -v numhit=${NUMHIT} -v maq=${MAQ} 'BEGIN{
  FS="\t"; OFS="\t";
  genepos = -1; geneposG = -1;
  cbpos = 0; umipos = 0; assignpos = 0;
  # for alignment with GN tag added
  cbposG = 0; umiposG = 0; assignposG = 0;
  correctCB = 0;
}{ 
  # for the first 1000 records, always check the position of CB or CR until CB has been identified
  # if the record has been assigned to a gene
  
  #### debug #######
  # print "before",$1,correctCB,cbpos,cbposG,assignpos,assignposG
  # if($1 != "J00151:455:HC57MBBXY:5:1117:12855:32068"){
  #   for(j = 1; j <= NF; j++){
  #     print j,$j
  #   }
  # }
  ##################
  
  if( ((correctCB == 0 && NR < 1000) || (cbposG == 0 || umiposG == 0)) && ($0 ~ /GN:Z:/) ){
    for(i = 12 ; i <= NF ; i++){
      curTag = substr($i,1,2)
      if( ( curTag == "CR" && correctCB == 0) || curTag == "CB"){
        cbposG = i
      }
      if( (curTag == "UR" && correctCB == 0) || curTag == "UB"){
        umiposG = i
      }
      if( (curTag == "XS" && correctCB == 0) || ($0 ~ /CB:Z:/ && curTag == "XS") ){
        assignposG = i
        geneposG = i + 2
      }
    }  
  }
  # if the record has not been assigned to a gene
  if( ((correctCB==0 && NR < 1000) || (cbpos ==0 || umipos == 0)) && !($0 ~/GN:Z:/) ){
	  for(i = 1 ; i <= NF ; i++){
	    curTag = substr($i,1,2)
      if( (curTag == "CR" && correctCB == 0) || curTag == "CB"){
        cbpos = i
        cbhead = curTag
      }
      if( (curTag == "UR" && correctCB == 0) || curTag == "UB"){
        umipos = i
        umihead = curTag
      }
      if((curTag == "XS" && currectCB == 0) || ($0 ~ /CB:Z:/ && curTag == "XS")){
        assignpos = i
        genepos = i + 2
      }
    }
  }
  
  # check if CB has been identified
  if(NR < 1000 && $0 ~ /CB:Z:/){
    correctCB = 1
  }
  # print "after",$1,correctCB,cbpos,cbposG,assignpos,assignposG
  # if CB has been corrected, then skip record if CB is not present
  if(correctCB == 1){
    if(substr($cbposG, 1, 2) != "CB" && substr($cbpos, 1, 2) != "CB"){
      # print cbposG,cbpos,$0
      next
    }
  }
  # extract information
  if ($assignposG == "XS:Z:Assigned"){
    # extract position of cb and umi
    nHit = substr($12,6) + 0
    if(nHit <= numhit){
      # cb, umi, gene; 
      cb = substr($cbposG, 6)
      umi = substr($umiposG, 6)
      gene = substr($geneposG, 6)
      count[cb"\t"umi"\t"gene]++
    }
  }
  if ($assignpos == "XS:Z:Assigned"){
    # extract position of cb and umi
    nHit = substr($12,6) + 0
    if(nHit <= numhit){
      # cb, umi, gene; 
      cb = substr($cbpos, 6)
      umi = substr($umipos, 6)
      gene = substr($genepos, 6)
      count[cb"\t"umi"\t"gene]++
    }
  }
}END{ 
  print cbhead"\t"umihead"\tgene\tcount"
  for(i in count){ 
    print i,count[i]; 
  }
}'
