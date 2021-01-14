
#' @title Algorithmic Framework for Representational Learning on Graphs
#'
#' @param data input data for edges consisting of at least two columns of data and if there are weights,it must be in the third column.
#' @param p return parameter.Default to 1.
#' @param q in-out parameter.Default to 1.
#' @param directed the network is directed or undirected.Default to undirected.
#' @param num_walks number of walks per node.Default to 10.
#' @param walk_length number of nodes in each walk.Default to 80.
#' @param dim embedding dimensions.Default to 128.
#'
#' @return embedding results for each node
#' @import  data.table
#' @import word2vec
#' @import rlist
#' @import vegan
#' @import vctrs
#' @importFrom dplyr %>%
#' @importFrom dplyr union
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph V
#' @importFrom igraph neighbors
#' @importFrom stats runif
#' @export
#' @examples
#' #Parameters can be customized as needed
#' data(gene_edges)
#' use_data<-gene_edges
#' emb<-node2vecR(use_data,p=2,q=1,num_walks=5,walk_length=5,dim=10)
#'
#'
#'
#'
#'

node2vecR<-function(data,p=NULL,q=NULL,directed=NULL,num_walks=NULL,walk_length=NULL,dim=NULL){

  if(is.null(directed)){
    directed="undirected"
  }
  if(is.null(p)){
    p<-1
  }
  if(is.null(q)){
    q<-1
  }
  if(is.null(num_walks)){
    num_walks<-10
  }
  if(is.null(walk_length)){
    walk_length<-80
  }
  if(is.null(dim)){
    dim<-128
  }
  edges<-data
  if(ncol(edges)==2){
    weights=rep(1,nrow(edges))
    edges1<-data.frame(edges,weight=weights)
  }else if(ncol(edges)==3)
  {
    weights=edges[,3]
    edges1<-data.frame(edges[,1:2],weight=weights)
  }
  names(edges1)<-c("no1","no2","weight")
  nodes<-union(edges[,1],edges[,2])
  G<- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

  get_nbr<-function(node)
  {
    node1<-as.character(node)
    node_nei<-V(G)$name[neighbors(G, node1, mode = "total")]
    node_nei<-data.frame(node_nei)
    names(node_nei)<-as.character(paste(node,"_nbr",sep=''))
    return(node_nei)
  }

  neigh<-list()
  for(i in 1:length(nodes)){
    node<-nodes[i]
    node_nbr<-get_nbr(node)
    neigh<-append(neigh,node_nbr)
  }


  node_probs<-function(node){
    weight<-NULL
    node_nei<-paste(node,"_nbr",sep='')
    if(directed=='undirected'){
      edge_weight<-edges1%>%select(weight)%>%filter((edges1$no1%in%neigh[[node_nei]]&edges1$no2==node)|
                                                      (edges1$no2%in%neigh[[node_nei]]&edges1$no1==node))
    }
    else{
      edge_weight<-edges1%>%select(weight)%>%filter(edges1$no2%in%neigh[[node_nei]]&edges1$no1==node)
    }
    edge_weight1<-as.numeric(t(edge_weight))
    total<-sum(edge_weight)
    probs<-edge_weight1/total
    probs1<-t(probs)
    probs2<-list(probs)
    names(probs2)<-paste(node,"_probs",sep='')
    return(probs2)
  }


  normalized_probs<-list()
  for(i in 1:length(nodes)){
    node<-nodes[i]
    b<-node_probs(node)
    normalized_probs<-append(normalized_probs,b)
  }
  print("The normalized probabilities finished.")


  alias_setup<-function(probs){
    K<-as.numeric(length(probs))
    q<-vector(mode="numeric",length=K)
    J<-vector(mode="numeric",length=K)
    smaller<-c()
    larger<-c()
    if(K>0){
      for(i in 1:K){
        q[i]<-(as.numeric(probs[i]))*K
        if(q[i]>1){
          larger<-append(larger,i)
        }else if(q[i]<1)
          smaller<-append(smaller,i)
      }
    }
    while(length(smaller)>0&&length(larger)>0){
      small<-smaller[length(smaller)]
      smaller<-smaller[-length(smaller)]
      large<-larger[length(larger)]
      larger<-larger[-length(larger)]
      J[small]<-large
      q[large]<-q[large]+q[small]-1
      if(q[large]<1){
        smaller<-append(smaller,large)
      }else{
        larger<-append(larger,large)
      }
    }
    alias_return<-rbind(q,J)
    return(alias_return)
  }


  alias_draw<-function(J,q){
    K=length(J)
    c<-runif(1,min=0,max=1)
    kk<-floor(c*K)
    if(kk==0){
      return(sample(1:length(q),1))
    }else{
      if(c<q[kk]){
        b<-kk
      }else{
        b<-J[kk]
      }
    }
    return(b)
  }

  node_alias_setup<-function(node){

    nbr<-get_nbr(node)
    probs<-node_probs(node)
    use_probs<-probs[[paste(node,"_probs",sep='')]]
    K<-as.numeric(length((use_probs)))
    result<-alias_setup(use_probs)
    return(result)
  }


  get_alias_1node_edge<-function(p,q,src,dst){
    no1<-NULL
    no2<-NULL
    weight<-NULL
    x<-NULL
    nei_dst<-unique(get_nbr(dst))  #dst's neighbor node
    nbr<-sort(nei_dst[,1])
    table<-c()
    for(i in 1:length(nbr)){
      dst_nbr<-nbr[i]
      total=0
      a<-edges1%>%select(no1,no2)%>%
        filter((edges1$no1==dst_nbr&edges1$no2==src))
      aim_weight<-edges1%>%select(weight)%>%
        filter((edges1$no1==dst&edges1$no2==dst_nbr))
      weight1<-as.numeric(unlist(aim_weight))
      if(length(weight1)!=0){
        if(dst_nbr==src){
          x<-weight1/p
        }else if(dst_nbr!=src&is.null(a)){
          x<-weight1/q
        }else if(is.null(a)==FALSE){
          x<-weight1
        }
      }
      y<-c(dst,dst_nbr,weight=x)
      if(length(y)==3){
        table<-rbind(table,y)
      }
    }
    table1<-as.data.frame(table)
    total<-0
    if(is.null(table1$weight)==FALSE){
      for(i in 1:nrow(table1)){
        total<-as.numeric(table1$weight[i])+total
      }
      for(i in 1:nrow(table1)){
        table1$weight[[i]]<-as.numeric(table1$weight[[i]])/total
      }
      n1<-as.character(table1[,1])
      n2<-as.character(table1[,2])
      n3<-as.character(table1$weight)
      K<-as.numeric(length(n3))
      n33<-alias_setup(n3)
      data.frame(n33)
    }else{
      n33<-NULL
    }
    return(n33)
  }

  alias_edges<-list()
  if(directed=="directed"){
    for(i in 1:nrow(edges1)){
      print(i)
      src<-edges1[i,]$no1
      dst<-edges1[i,]$no2
      if(length(get_nbr(dst))>0){
        y<-get_alias_1node_edge(p,q,src,dst)
        y1<-list(y)
        names(y1)<-paste(src,"-",dst,sep='')
      }else break
      alias_edges<-append(alias_edges,y1)
    }
  }else if(directed=="undirected"){
    for(i in 1:nrow(edges1)){
      print(i)
      src<-edges1[i,]$no1
      dst<-edges1[i,]$no2
      if(length(get_nbr(dst))>0){
        y<-get_alias_1node_edge(p,q,src,dst)
        y1<-list(y)
        z<-get_alias_1node_edge(p,q,dst,src)
        z1<-list(z)
        names(y1)<-paste(src,"-",dst,sep='')
        names(z1)<-paste(dst,"-",src,sep='')
      }else break
      alias_edges<-append(alias_edges,y1)
      alias_edges<-append(alias_edges,z1)
    }
  }



  alias_nodes<-list()
  for(i in 1:length(nodes)){

    node<-nodes[i]
    alias_node<-node_alias_setup(node)
    #r<-data.frame(node_name=node,probs=re1,alias=re2)
    alias_node1<-list(name=alias_node)
    names(alias_node1)<-paste(node,"-alias",sep='')
    alias_nodes<-append(alias_node1,alias_nodes)
  }


  node2vec_walk<-function(walk_length,start_node){
    walk<-c()
    walk<-append(start_node,walk)
    while(length(walk)<walk_length){
      cur<-walk[length(walk)]
      cur_nbrs<-sort(get_nbr(cur)[,1])

      if(length(cur_nbrs)>0){
        if(length(walk)==1){
          cur_alias<-paste(cur,"-alias",sep='')
          J<-(alias_nodes[[cur_alias]])[2,]
          q<-(alias_nodes[[cur_alias]])[1,]
          n<-alias_draw(J,q)
          if(n<=length(cur_nbrs)){
            walk1<-cur_nbrs[n]
            walk<-append(walk,walk1)
          }else{
            walk<-walk
          }
        }
        else if(length(walk)>1){
          prev<-walk[length(walk)-1]
          aim_edge<-alias_edges[[paste(prev,"-",cur,sep='')]]
          J<-aim_edge[2,]
          q<-aim_edge[1,]
          n<-alias_draw(J,q)
          if(n<=length(cur_nbrs)){
            next_node<-cur_nbrs[alias_draw(J,q)]
            walk<-append(walk,next_node)
          }else{
            walk<-walk
          }
        }
      }

    }
    return(walk)
  }

  simulate_walk<-function(num_walks,walk_length){
    print("walk interation:")
    final_walk<-list()
    for(walk_iter in 1:num_walks){
      output<-paste(as.numeric(walk_iter),"/",as.numeric(num_walks),sep='')
      print(output)
      use_nodes<-sample(nodes,length(nodes))
      for(i in 1:length(use_nodes)){
        walk<-node2vec_walk(walk_length,start_node=use_nodes[i])
        final_walk<-append(walk,final_walk)
      }
    }
    return(final_walk)
  }

  walks<-simulate_walk(num_walks,walk_length)
  walk_use<-as.character(walks)
  dims<-dim
  model<-word2vec(walk_use,type="skip-gram",dim=dims)
  emb<-as.matrix(model)
  emb_final<-emb[!rownames(emb) %in% c("</s>"),]
  return(emb_final)
}




