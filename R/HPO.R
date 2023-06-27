hpo_github_base_url <-'https://github.com/obophenotype/human-phenotype-ontology/'
inheritance_term_id <- 'HP:0000005'

## TODO: move these to options in exported functions(e.g. secure = TRUE), use as httr::with_config
#'@export
insecure <- function() httr::set_config(httr::config(ssl_verifypeer = 0L))
#'@export
secure <- function() httr::set_config(httr::config(ssl_verifypeer = 1L))

latest_hpo_release <-function() {
  (function() {
    # attempt to get latest hpo release from github, but otherwise use cached version in case server is down
    tag_url <- 'https://api.github.com/repos/obophenotype/human-phenotype-ontology/tags'
    
    tryCatch(
      content(retry(tag_url, verb = 'GET')) %>% 
        map_chr('name') %>% 
        keep(str_detect, 'v[0-9]{4}-[0-9]{2}-[0-9]{2}') %>% 
        sort(decreasing = TRUE) %>% 
        first() %>% 
        (function(x) { assert_that(is_scalar_character(x) & !is.na(x)); x }),
      error = function(e) {
        hpo_files <- list.files(get_cache_dir(), pattern = '^hpo\\..*\\v[0-9]{4}-[0-9]{2}-[0-9]{2}\\.rds$')
        if (length(hpo_files)) {
          ver <- 
            str_extract(hpo_files, 'v[0-9]{4}-[0-9]{2}-[0-9]{2}(?=\\.rds)') %>% 
            sort(decreasing = TRUE) %>% 
            first()
          warning("Coudn't retrieve latest HPO release from: ", tag_url, '. ',
                  "Using cached version ", ver, '.')
          ver
        } else {
          stop('could not get latest hpo release')
        }
      })
  }) %>% cache('latest_hpo_release')
} 

# this contains childmost term for a phenotype heirachy
# use for annotating genes
get_genes_to_phenotype <- function() 
{
  hpo_release <- latest_hpo_release()
  
  url <- str_c(hpo_github_base_url, 'releases/download/', hpo_release, '/genes_to_phenotype.txt')
  
  (function()
    retry('GET', url) %>% 
      content(as = 'raw') %>% 
      rawConnection() %>% 
      read_tsv(col_names = c('entrez_gene_id',
                             'entrez_gene_symbol',
                             'hpo_term_id',
                             'hpo_term_name', 
                             'frequency', 
                             'disease_id'),
               skip = 1,
               col_types = 'iccccc')) %>% 
    cache('hpo.genes_to_phenotype',
          ver = hpo_release,
          disk = TRUE)
}

# this contains all hpo terms
# use for creating gene lists from hpo terms
get_phenotype_to_genes <- function() 
{
  hpo_release <- latest_hpo_release()
  
  url <- str_c(hpo_github_base_url, 'releases/download/', hpo_release, '/phenotype_to_genes.txt')
  
  (function()
    retry('GET', url) %>% 
      content(as = 'raw') %>% 
      rawConnection() %>% 
      read_tsv(col_names = c('hpo_term_id',
                             'hpo_term_name',
                             'entrez_gene_id', 
                             'entrez_gene_symbol',
                             'disease_id'),
               col_types = 'ccicc',
               skip = 1)) %>% 
    cache('hpo.phenotype_to_genes',
          ver = hpo_release,
          disk = TRUE)
}

get_omim_gene_map <- function() 
{
  (function()
    get_genes_to_phenotype() %>% 
     select(entrez_gene_id, entrez_gene_symbol, disease_id, hpo_term_id) %>% 
     filter(str_starts(disease_id, 'OMIM:')) %>% 
     left_join(select(get_inheritance_terms(), hpo_term_id) %>% 
                 mutate(is_inh = TRUE),
               by = "hpo_term_id") %>% 
     mutate(is_inh = replace_na(is_inh, FALSE)) %>% 
     group_by(entrez_gene_id, entrez_gene_symbol, disease_id) %>% 
     summarise(inheritance = list(hpo_term_id[is_inh]),
               .groups = 'drop') %>% 
     mutate(inheritance = simplify_inheritance(inheritance)) %>% 
     mutate(hgnc_symbol = coalesce(hgnc_entrez2sym(entrez_gene_id),
                                   hgnc_sym2sym(entrez_gene_symbol),
                                   entrez_gene_symbol)) %>% 
     select(symbol = hgnc_symbol, disease_id, inheritance)) %>% 
    cache('hpo_omim_gene_map',
          disk = TRUE,
          ver = latest_hpo_release())
}

get_term_names <- function() 
{
  (function()
    get_phenotype_to_genes() %>% 
     select(hpo_term_id, hpo_term_name) %>% 
     distinct() %>% 
     arrange(hpo_term_id)) %>% 
    cache('term_names')
}

#'@importFrom httr GET accept_json content
#'@importFrom stringr str_c
hpo_api_get <- function(extension,
                        base_url = "https://hpo.jax.org/api/hpo/")
{
  url <- str_c(base_url, extension)
  
  response <- retry('GET', url, accept_json())
  
  if (response$status_code == 200) {
    return(content(response))
  }
  return(NULL)
}

#'@importFrom dplyr tibble as_tibble bind_rows "%>%"
#'@importFrom memoise memoise 
#'@importFrom purrr map_df
#'@export
get_hpo_term <- function(hpo_id)
{
  assert_that(is.character(hpo_id),
              all(str_detect(hpo_id, 'HP:\\d{7}'), na.rm = TRUE))
  
  mapper <- cavalier_cache$hpo_get_term_mapper
  
  if (is.null(mapper)) {
    mapper <- memoise(
      function(hpo_id) {
        if (!is.na(hpo_id)) {
          result <- hpo_api_get(str_c('term/', hpo_id))
          if (!is.null(result)) {
            data <-
              result$details %>%
              map(~ `if`(is.list(.), list(unlist(.)), .)) %>% 
              as_tibble() %>% 
              mutate(parents = list(bind_rows(result$relations$parents)),
                     children = list(bind_rows(result$relations$children))) %>% 
              rename(ontologyId = id) %>% 
              select(ontologyId, everything())
            return(data)
          } 
        }
        return(tibble(ontologyId = hpo_id))
      })
    cavalier_cache$hpo_get_term_mapper <- mapper
  }
  
  map_df(hpo_id, mapper)
}

#'@importFrom rlang is_integerish 
#'@export
get_hpo_term_genes <- function(hpo_id)
{
  assert_that(is.character(hpo_id),
              all(str_detect(hpo_id, 'HP:\\d{7}'), na.rm = TRUE))
  
  mapper <- cavalier_cache$get_hpo_term_genes_mapper
  
  if (is.null(mapper)) {
    mapper <- memoise(
      function(hpo_id) {
        if (!is.na(hpo_id)) {
          result <- hpo_api_get(str_c('term/', hpo_id, '/genes?max=-1'))
          if (!is.null(result)) {
            data <-
              tibble(ontologyId = hpo_id,
                     genes = list(
                       map_df(result$genes,
                              ~ tibble(entrezGeneId = .$entrezGeneId,
                                       entrezGeneSymbol = .$entrezGeneSymbol,
                                       dbDiseases = list(bind_rows(.$dbDiseases))))
                     ))
            return(data)
          } 
        }
        return(tibble(ontologyId = hpo_id))
      })
    cavalier_cache$get_hpo_term_genes_mapper <- mapper
  }
  
  map_df(hpo_id, mapper)
}

#'@importFrom rlang is_integerish 
#'@export
get_hpo_gene <- function(entrez_id)
{
  assert_that(is_integerish(entrez_id),
              all(entrez_id > 0, na.rm = TRUE))
  
  mapper <- cavalier_cache$get_hpo_gene_mapper
  
  if (is.null(mapper)) {
    mapper <- memoise(
      function(entrez_id) {
        if (!is.na(entrez_id)) {
          result <- hpo_api_get(str_c('gene/', entrez_id))
          if (!is.null(result)) {
            data <- 
              as_tibble(result$gene) %>% 
              mutate(termAssoc = list(bind_rows(result$termAssoc)),
                     diseaseAssoc = list(bind_rows(result$diseaseAssoc)))
            return(data)
          } 
        }
        return(tibble(entrezGeneId = entrez_id))
      })
    cavalier_cache$get_hpo_gene_mapper <- mapper
  }
  
  # map_df_prog(entrez_id, mapper)
  map_df(entrez_id, mapper)
}

#'@importFrom purrr map_int
get_term_descendants <- function(term_id) 
{
  recursion <- cavalier_cache$get_term_descendants_recursion
  
  if (is.null(recursion)) {
    recursion <- memoise(function(init) {
      init %>%
        filter(childrenCount > 0) %>%
        pull(ontologyId) %>%
        map_df(function(id) {
          get_hpo_term(id) %>% 
            select(children) %>% 
            unnest(children) %>% 
            mutate(parentId = !!id) %>% 
            select(ontologyId, name, parentId, childrenCount) %>% 
            recursion()
        }) %>% 
        bind_rows(init, .)
    })
    cavalier_cache$get_term_descendants_recursion <- recursion
  }
  
  get_hpo_term(term_id) %>% 
    mutate(parentId = NA_character_,
           childrenCount = map_int(children, nrow)) %>% 
    select(ontologyId, name, parentId, childrenCount) %>% 
    recursion()
}

get_inheritance_terms <- function() 
{
  (function()
    get_term_descendants(inheritance_term_id) %>% 
     select(hpo_term_id = ontologyId,
            hpo_term_name = name,
            parent_term_id = parentId,
            num_children = childrenCount)) %>% 
    cache('inheritance_terms')
}

#' #'@export
#' get_hpo_api_disease <- function(disease_id)
#' {
#'   assert_that(
#'     is.character(disease_id),
#'     all(str_detect(disease_id, '^(OMIM)|(ORPHA):\\d+$'), na.rm = TRUE))
#'   
#'   mapper <- cavalier_cache$get_hpo_disease_mapper
#'   
#'   if (is.null(mapper)) {
#'     mapper <- memoise(
#'       function(disease_id) {
#'         if (!is.na(disease_id)) {
#'           result <- hpo_api_get(str_c('disease/', disease_id))
#'           if (!is.null(result)) {
#'             data <-
#'               as_tibble(result$disease) %>% 
#'               mutate(
#'                 catTermsMap = list(map_df(result$catTermsMap, function(data) {
#'                   tibble(catLabel = data$catLabel,
#'                          terms = list(bind_rows(data$terms))) %>% 
#'                     unnest(terms)
#'                 })),
#'                 geneAssoc = list(bind_rows(result$geneAssoc)))
#'             return(data)
#'           } 
#'         }
#'         return(tibble(diseaseId = disease_id,
#'                       diseaseName = NA_character_))
#'       })
#'     cavalier_cache$get_hpo_disease_mapper <- mapper
#'   }
#'   
#'   # map_df_prog(disease_id, mapper)
#'   map_df(disease_id, mapper)
#' }

# pull disease names from HPO API and store locally
#' @export
disease_names <- function(ids)
{
  hpo_release <- latest_hpo_release()
  
  dndb_0 <-
    (function() 
      tibble(disease_id   = character(), 
             disease_name = character())) %>% 
    cache('hpo_disease_names',
          disk = TRUE,
          ver = hpo_release)
  
  dndb <- 
    tibble(disease_id = ids) %>% 
    filter(str_detect(disease_id, '(OMIM)|(ORPHA):[0-9]+')) %>% 
    anti_join(dndb_0,  by = "disease_id") %>% 
    mutate(disease_name = map_chr(disease_id, function(x) {
      result <- hpo_api_get(str_c('disease/', x))$disease$diseaseName
      `if`(is.null(result), NA_character_, result)
    })) %>% 
    bind_rows(dndb_0) %>% 
    arrange_all()
  
  if (nrow(dndb) > nrow(dndb_0)) {
    (function() dndb) %>% 
      cache('hpo_disease_names',
            disk = TRUE,
            ver = hpo_release,
            overwrite = TRUE)
  }
  
  with(dndb, disease_name[match(ids, disease_id)])
}

#' @export
#' @importFrom dplyr inner_join anti_join add_row
get_hpo_gene_list <- function(hpo_id, prefer_omim = TRUE) {
  
  assert_that(is_scalar_character(hpo_id),
              !is.na(hpo_id),
              str_detect(hpo_id, 'HP:\\d{7}'))
  
  term_name <- 
    with(get_term_names(),
         hpo_term_name[match(hpo_id, hpo_term_id)])
  
  get_phenotype_to_genes() %>% 
    filter(hpo_term_id == hpo_id) %>% 
    select(entrez_gene_id, entrez_gene_symbol, disease_id) %>%
    group_by(entrez_gene_id) %>% 
    filter(!prefer_omim | str_starts(disease_id, 'OMIM') | !any(str_starts(disease_id, 'OMIM'))) %>% 
    ungroup() %>% 
    (function(gd) {
      gd %>% 
        inner_join(get_genes_to_phenotype() %>% 
                     select(hpo_term_id, entrez_gene_id, disease_id),
                   by = c("entrez_gene_id", "disease_id")) %>% 
        inner_join(select(get_inheritance_terms(), hpo_term_id),
                   by = "hpo_term_id") %>% 
        rename(inheritance = hpo_term_id) %>% 
        chop(inheritance) %>% 
        complete(gd) %>% 
        mutate(inheritance = simplify_inheritance(inheritance))
    }) %>% 
    mutate(gene = coalesce(hgnc_entrez2sym(entrez_gene_id),
                           hgnc_sym2sym(entrez_gene_symbol),
                           entrez_gene_symbol)) %>% 
    select(gene, disease_id,  inheritance) %>%
    arrange(gene, disease_id) %>% 
    mutate(., version = str_c('build ', latest_hpo_build_num())) %>% 
    mutate(list_id = hpo_id,
           list_name = term_name) %>% 
    select(list_id, list_name, everything())
}

hpo_simple_inheritnace <-
  tibble(name = character(), ontologyId = list()) %>% 
  add_row(name = 'AR', ontologyId = list('HP:0000007')) %>% 
  add_row(name = 'AD', ontologyId = list('HP:0000006')) %>% 
  add_row(name = 'AD', ontologyId = list('HP:0025352')) %>% 
  add_row(name = 'XL', ontologyId = list('HP:0001417')) %>% 
  add_row(name = 'XL', ontologyId = list('HP:0001419')) %>% 
  add_row(name = 'XL', ontologyId = list('HP:0001423')) %>% 
  add_row(name = 'XL', ontologyId = list(c('HP:0001417', 'HP:0001419'))) %>% 
  add_row(name = 'XL', ontologyId = list(c('HP:0001417', 'HP:0001423'))) %>% 
  add_row(name = 'XL', ontologyId = list(c('HP:0001419', 'HP:0001423'))) %>% 
  add_row(name = 'MT', ontologyId = list('HP:0001427')) %>% 
  add_row(name = 'MT', ontologyId = list(c('HP:0000007', 'HP:0001427'))) %>% 
  add_row(name = 'AR/AD', ontologyId = list(c('HP:0000006', 'HP:0000007')))

hpo_inheritnace_ignore <- c(
  'HP:0001425', # Heterogeneous 
  'HP:0003745', # Sporadic,
  'HP:0003743', # Genetic anticipation
  'HP:0001428', # Somatic mutation
  'HP:0001426', # Multifactorial inheritance
  'HP:0001452', # Autosomal dominant contiguous gene syndrome
  'HP:0001466', # Contiguous gene syndrome
  'HP:0001442', # Somatic mosaicism
  'HP:0010984', # Digenic inheritance
  'HP:0010982', # Polygenic inheritance
  'HP:0003744'  # Genetic anticipation with paternal anticipation bias
)

simplify_inheritance <- function(inh)
{
  tibble(ontologyId = as.list(inh)) %>% 
    mutate(ontologyId = map(ontologyId, ~ sort(setdiff(., hpo_inheritnace_ignore)))) %>%
    left_join(hpo_simple_inheritnace, by = "ontologyId") %>% 
    mutate(name = replace_na(name, 'OTHER')) %>%
    pull(name)
}
