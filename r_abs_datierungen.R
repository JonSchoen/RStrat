# --- --- --- --- --- FUNKTION absolute_datierungen --- --- --- --- --- --- ----

#' Adds absolute datings
#'
#' @export
#' @param absolute_data_list PLATZHALTERPLATZHALTER
#' @param liste_stellen PLATZHALTERPLATZHALTER
#' @param matrix_gleich PLATZHALTERPLATZHALTER
#' @param matrix_ueber_unter PLATZHALTERPLATZHALTER
#' @return PLATZHALTER PLATZHALTERPLATZHALTERPLATZHALTER
#' @examples
#' absolute_chronology <- absolute_datierungen(absolute_data_list, liste_stellen, matrix_gleich, matrix_ueber_unter)




absolute_datierungen <- function(absolute_data_list, liste_stellen, matrix_gleich, matrix_ueber_unter)
  # liste_stellen ist redundant, könnnte noch wegrationalisiert werden

{
  # --- --- --- 1.1. Aufbereitung der Daten --- --- ----

  absolute_data_list[,1] <- as.character(absolute_data_list[,1])
  absolute_data_list[,2] <- as.double(absolute_data_list[,2])


  # Abgleichen, ob alle Befunde aus der Liste der absoluten Datierungen auch in den Matritzen vorhanden sind
  break_req = FALSE

  for(fortl_nr in 1:nrow(absolute_data_list))
  {
    stelle <- absolute_data_list[fortl_nr,1]

    if (stelle %in% liste_stellen == FALSE)
    {
      conflict <- paste("Feature ", absolute_data_list[fortl_nr,1], " in your list of absolute datiengs does not appear in your stratification and must be included.")
      report_of_conflict <- list(break_req = TRUE, conflict = conflict)

      analysis_functions_objects()

      return(report_of_conflict)
    }# Ende if-loop
  }# ende for-loop

  #Erstellen einer Liste mit allen Stellen und deren absoluten Datierungen
  absolute_data_list_gesamt <- liste_alle_befunde_mit_dat(liste_stellen, absolute_data_list, matrix_gleich)

  # --- --- --- 1.2. Erstellen der Matritzen --- --- ----
  #Zusammenfuehren der Matrix und der Liste absoluter Datierungen
  #Anh?ngen einer Zeile/Spalte an die Matritzen f?r Datierung
  matrix_ueber_unter_daten <- cbind(date=NA,matrix_ueber_unter)
  matrix_ueber_unter_daten <- rbind(date=NA,matrix_ueber_unter_daten)

  # Einf?gen der absoluten Datierungen in die Matrix
  matrix_ueber_unter_daten[1,2:ncol(matrix_ueber_unter_daten)] <- absolute_data_list_gesamt[,2]
  matrix_ueber_unter_daten[2:nrow(matrix_ueber_unter_daten),1] <- absolute_data_list_gesamt[,2]


  # --- --- --- 2. Loeschen redundanter Daten --- --- ----
  # Loeschen der absoluten Datierungen, die irrelevant sind, da ein j?ngeres Datum in einer unteren Schicht vorhanden ist

  matrix_ueber_unter_daten2 <- matrix_ueber_unter_daten
  #Loeschen aller NA's
  NA_index <- which(is.na(matrix_ueber_unter_daten2[1,2:ncol(matrix_ueber_unter_daten2)]))+1   #gibt Vektor mit Positionen aller NA's in den Spalten.
  #2:ncol, da die erste Spalte die Datierungen enth?lt.
  #"+1", da wegen 2:ncol die zweite Spalte als erste gez?hlte wird etc.
  matrix_ueber_unter_daten2 <- matrix_ueber_unter_daten2[,-NA_index]
  matrix_ueber_unter_daten2 <- t(matrix_ueber_unter_daten2)
  NA_index <- which(is.na(matrix_ueber_unter_daten2[1,2:ncol(matrix_ueber_unter_daten2)]))+1   #gibt Vektor mit Positionen aller NA's in den Spalten.
  #2:ncol, da die erste Spalte die Datierungen enth?lt.
  #"+1", da wegen 2:ncol die zweite Spalte als erste gez?hlte wird etc.
  matrix_ueber_unter_daten2 <- matrix_ueber_unter_daten2[,-NA_index]
  matrix_ueber_unter_daten2 <- t(matrix_ueber_unter_daten2)


  # erzeugt eine Matrix in der die redundanten Daten aufgelistet werden
  redundante_daten <- matrix(c("Das Datum", "aus dem Befund","liegt stratigrafisch ?ber dem Datum","aus dem Befund"), nrow = 1, ncol = 4)

  # stellt redundante Daten fest und traegt diese in die Matrix ein
  for(row in 2:nrow(matrix_ueber_unter_daten2))
  {
    for(col in 2:ncol(matrix_ueber_unter_daten2))
    {
      if(matrix_ueber_unter_daten2[row,col] > 0 && matrix_ueber_unter_daten2[row,1] < matrix_ueber_unter_daten2[1,col])
      {
        absolute_data_list_gesamt[which(absolute_data_list_gesamt[,2] == matrix_ueber_unter_daten2[row,1]),2] <- NA # [which(...),2] gibt die Koordinate mit dem Datum an, das auf NA gesetzt werden soll

        redudata_temp <- c(matrix_ueber_unter_daten2[row,1], rownames(matrix_ueber_unter_daten2)[row], matrix_ueber_unter_daten2[1,col], colnames(matrix_ueber_unter_daten2)[col] )
        redundante_daten <- rbind(redundante_daten, redudata_temp)# bei jedem Durchlauf wird eine neue Zeile angehaengt

      }
    }
  }

#  write.table(redundante_daten, file="C:\\Users\\Johnny B.\ Goode\\Desktop\\R\\2017_02_27\\Output\\redundante_daten.csv", row.names=FALSE, col.names=FALSE, sep = " , ")
  print(paste("Redundante Datierungen fallen gelassen"))

  #Erneutes bauen der Matrix, ohne Datierungen die aufgrund ihrer stratigraphischen Einordnung hinf?llig sind
  matrix_ueber_unter_daten[1,2:ncol(matrix_ueber_unter_daten)] <- absolute_data_list_gesamt[,2]
  matrix_ueber_unter_daten[2:nrow(matrix_ueber_unter_daten),1] <- absolute_data_list_gesamt[,2]


  # --- --- --- 3.1. sorteiren der Matrix nach Datierung --- --- ----
  matrix_nach_daten <- matrix_ueber_unter_daten[order(matrix_ueber_unter_daten[,1], na.last = FALSE),]
  matrix_nach_daten <- t(matrix_nach_daten)

  matrix_nach_daten <- matrix_nach_daten[order(matrix_nach_daten[,1], na.last = FALSE),]
  matrix_nach_daten <- t(matrix_nach_daten)


  # --- --- --- 3.2. Splitten und Kuerzen --- --- ----
  # Splitten in 2 Matritzen gleicher Leserichtung
  matrix_ueber <- matrix_nach_daten
  matrix_unter <- t(matrix_nach_daten)



  # Loeschen absolut datierter Zeilen und nicht datierter Spalten
  matrix_ueber <- kuerzen_matrix(matrix_ueber)
  matrix_unter <- kuerzen_matrix(matrix_unter)



  # --- --- --- 4. Erstellen einer Tabelle mit absoluten Datierungen --- --- ----
  tabelle_dat <- erstelle_tabelle(matrix_ueber, matrix_unter, absolute_data_list_gesamt)
  absolute_date <- list(break_req = break_req, tabelle_dat = tabelle_dat)

  analysis_functions_objects()

  return(absolute_date)
}
