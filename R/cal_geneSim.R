##' @title Calculate the functional similarities between two genes based on Gene Ontology.
##' @description Returns the pairwise similarities between two genes. Different calculation geneFun are provided as implemented in GOSim R package.
##' @param gene1 The first gene, such as "KY.Chr1.1".
##' @param gene2 The second gene, such as "KY.Chr1.10".
##' @param geneID2GO A list contains gene IDs and GO terms.
##' @param simFun The calculation method between pairwise GO terms.
##' @param geneFun The calculation method between pairwise genes.
##' @param category Gene Ontology categories. Cellular component as CC; biological process as BP; molecular function as MF.
##' @param return This function will stop by default when error happens. If TRUE, then it will return NA.
##' @return A value of functional similarity.
##' @examples
##' cal_geneSim("KY.Chr1.1","KY.Chr1.1000")
##' @author Chen Zaohuang
##' @import topGO
##' @import GOSim
##' @export

cal_geneSim<-function(gene1,gene2,geneID2GO=KYgeneID2GO,simFun=c("Resnik","Lin","JiangConrath","relevance"),geneFun=c("mean","max","simAvg","simMax","simMin"),category="BP",return=F){

  simFun<-match.arg(simFun,choices = c("Resnik","Lin","JiangConrath","relevance"))
  geneFun<-match.arg(geneFun,choices = c("mean","max","simAvg","simMax","simMin"))

    if(all(c(gene1,gene2) %in% names(geneID2GO)) ){

      gene1.GO<-unlist(geneID2GO[names(geneID2GO) == gene1],use.names = F)
      gene2.GO<-unlist(geneID2GO[names(geneID2GO)==gene2],use.names = F)

      if( category=="BP"){
        BPterms <- ls(GOBPTerm)
        data("ICsBPhumanall")
        gene1.GO <- gene1.GO[gene1.GO %in% BPterms & gene1.GO %in% names(IC)]
        gene2.GO <- gene2.GO[gene2.GO %in% BPterms & gene2.GO %in% names(IC)]


      } else if(category=="MF"){
        MFterms<-ls(GOMFTerm)
        data("ICsMFhumanall")
        gene1.GO <- gene1.GO[gene1.GO %in% MFterms & gene1.GO %in% names(IC)]
        gene2.GO <- gene2.GO[gene2.GO %in% MFterms & gene2.GO %in% names(IC)]

      } else {
        CCterms<-ls(GOCCTerm)
        data("ICsCChumanall")
        gene1.GO <- gene1.GO[gene1.GO %in% CCterms & gene1.GO %in% names(IC)]
        gene2.GO <- gene2.GO[gene2.GO %in% CCterms & gene2.GO %in% names(IC)]

      }

      unique(gene1.GO)->gene1.GO
      unique(gene2.GO)->gene2.GO

      if ( length(gene1.GO)>0 & length(gene2.GO) >0 ){

      }  else { return(NA)}

      df<-matrix(ncol = length(gene2.GO),nrow = length(gene1.GO))
      colnames(df)<-gene2.GO
      rownames(df)<-gene1.GO

      for (i in 1:length(gene1.GO)){
        for(j in 1:length(gene2.GO)){
          gosim<-getTermSim(c(gene1.GO[i],gene2.GO[j]),method = simFun,verbose = FALSE)
          df[i,j]<-gosim[upper.tri(gosim)]
        }
      }
      if(geneFun=="mean"){
        return(mean(df))
      } else if (geneFun=="max"){
        return(max(df))
      } else if ( geneFun=="SimAvg") {
        rowMax=mean(apply(df, 1, max))
        colMax=mean(apply(df, 2, max))
        return( 0.5 *(rowMax+colMax))
      } else if (geneFun=="SimMax"){
        rowMax=mean(apply(df, 1, max))
        colMax=mean(apply(df, 2, max))
        return(max(rowMax,colMax))
      }else if (geneFun=="SimMin"){
        rowMax=min(apply(df, 1, max))
        colMax=min(apply(df, 2, max))
        return(min(rowMax,colMax))
      }
    } else if (!(all(c(gene1,gene2) %in% names(geneID2GO))) & return==T ){
      return(NA)
    } else {stop ("Not all genes have GO terms") }



}



##' @title Get GO terms
##' @description Returns GO terms based the gene ID provided.
##' @param geneID Gene ID such as "KY.Chr1.1".
##' @param geneID2GO A list contains gene IDs and GO terms.
##' @param remove Keep gene ID or not. If True, then remove the geneID.
##' @return A vector of GO terms.
##' @examples
##' test<-get_GO("KY.Chr1.1",geneID2GO=KYgeneID2GO)
##' @author Chen Zaohuang
##' @export

get_GO<-function(geneID,geneID2GO=KYgeneID2GO,remove=T){

  geneGO <-geneID2GO[names(geneID2GO) %in% geneID]
  if(remove==T){
    geneGO<-unlist(unlist(geneGO,use.names = F))
    return(geneGO)
  } else {return(geneGO)}


}
