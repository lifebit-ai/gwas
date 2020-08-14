#' Title
#'
#' @param df 
#' @param nDigits_after_decimal 
#' @param table_caption 
#' @param escape 
#' @param font_family 
#'
#' @return
#' @export
#'
#' @examples
DTable <- function(df = df,
                   nDigits_after_decimal = 1,
                   table_caption = ""  ,
                   escape = FALSE,
                   font_family = "tahoma"){ # absolut path to file
  
  
  # > READING FILE INTO DF; DF TRANSFORMATIONS:
  # Read file
  
  # > Customize interactive DT::datatable
  DT::datatable(df, 
                rownames   = FALSE,
                escape     = FALSE,
                filter     = "bottom",
                caption    = table_caption,
                extensions = list("ColReorder"   = NULL,
                                  "Buttons"      = NULL,
                                  "FixedColumns" = list(leftColumns=1),
                                  "RowGroup"     = list(dataSrc="")),
                
                # OPTIONS:
                options = list(
                  
                  # Does not allow columnful dataframes go rogue and tucks them in to fit page width
                  scrollX = TRUE,
                  scrollCollapse = TRUE, 
                  
                  # Defines all capabilities
                  dom        = 'Brltpi', 
                  
                  autoWidth  = TRUE,
                  ColReorder = TRUE,
  
                  #   columnDefs = list(list(targets = length(colnames(df)), visible = TRUE)))),
                  
                  lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),
                  
                  buttons    = list('copy','print', 
                                    list(extend  = 'collection',
                                         buttons = c('csv', 'excel', 'pdf'),
                                         text    = 'Save'),
                                    I('colvis')),
                  # Black header container for colnames
                  initComplete = JS(
                    "function(settings, json) {",
                    "$(this.api().table().header()).css({'color': '#fff','font-family': 'tahoma', 'background-color': '#71879d'});",
                    "$(this.api().table().body()).css({'color': '#71879d','font-family': 'tahoma',   'text-align' : 'center'});",
                    "$(this.api().table().footer()).css({'color': '#fff','font-family': 'tahoma'});",
                    "$(this.api().table().container()).css({'color': '#fff','font-family': 'tahoma', 'outline-color' : '#71879d' });",
                    "$(this.api().table().node()).css({'color': '#fff','font-family': 'tahoma'});",
                    "}") )) %>% 
    
    # Change fontsize of cell values
    formatStyle(columns    = seq_along(colnames(df)), 
                fontSize   = "85%",
                fontFamily = "tahoma")%>%
    
    # Round values to 4 decimals
    formatRound(colnames(df), nDigits_after_decimal)  -> fancyDatatable
  
  return(fancyDatatable)  
}
