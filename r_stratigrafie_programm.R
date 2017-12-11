#' Stratifying archaeological data
#'
#' \code{RStratigraphy} returns different tables in which logical relations between a stratified data are completed, as well as their absolute datings. In addition, logical mistakes in the dataset are revealed.
#'
#' @export
#' @param strat_list A table of three columns and n rows, containing all known relations between the features. The first row has to name the columns. The first and the third columns contains the names of the features, where as the second column containes the relation in form of "above", "under" or "equal".
#' @param absolute_data_list A table of two columns and n rows, containing all known absolute datings of dated features. The first row has to name the columns. The first column contains the name of the feature, the second one the datings as numeric number.
#' @param fehler_loeschen A function, that automaticly erases logical mistakes in your stratigraphy if set on TRUE. WARNING! This function does not establish scientifically correct data!
#' @return A list of tables: under-above-relation of the features = tab_under_above, equal features = tab_equal, absolute datings of features = absolute_chronology.
#' @examples
#' RStratigraphy <- function(strat_list, absolute_data_list, fehler_loeschen)


RStratigraphy <- function(strat_list, absolute_data_list, fehler_loeschen)
{

  # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  # --- --- --- --- --- Formatierung der Daten fuer Weiterverarbeitung --- evtl. in separate Funktion--- --- ----

  # --- --- --- --- --- --- --- --- --- --- --- ---

  # Liste in Typ "character" umwandeln, da numerische und nicht numerische Befundbezeichnungen sonst nicht gleich behandelt werden:
  strat_list[] <- lapply(strat_list, as.character)

  # Erstellen einzelner Listen nach Lageverhaeltnis der Befunde untereinander:
  # Die Temp-Tabellen werden von der Fehlerkontrolle verwendet!
  liste_ueber_temp <- subset(strat_list, strat_list[,2] == "ueber")
  liste_unter_temp <- subset(strat_list, strat_list[,2] == "unter")
  liste_gleich_temp <- subset(strat_list, strat_list[,2] == "gleich")

  # Listen vervollstaendigen und fuer die Befuellung der Matrix vorbereiten
  # Unter-Liste in Ueber-Liste uebersetzen und mit bestehender ueber-Liste zusammenfuehren:
  liste_ueber_temp <- ergaenzen_listen(liste_ueber_temp, liste_unter_temp)
  # in Gleich-Liste werden umgekehrte Ausdruecke eingetragen (3 = 5 und 5 = 3):
  liste_gleich_temp <- ergaenzen_listen(liste_gleich_temp)

  # Spalte mit ueber/unter bzw. gleich loeschen:
  liste_ueber <- liste_ueber_temp[,c(1,3)]
  liste_gleich <- liste_gleich_temp[,c(1,3)]

  # Liste mit allen beteiligten Stellen erstellen
  # Liste Stellen wird aus Spalte links und Spalte rechts erstellt:
  liste_stellen <- c(array(strat_list[,1]),array(strat_list[,3]))

  # Sortieren der Stellen:
  liste_stellen <- liste_stellen[order(liste_stellen)]

  # Loeschen doppelt genannte Stellen.
  liste_stellen <- unique(liste_stellen)

  # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  # --- --- --- --- --- Erstellen der Matrizen --- evtl. in separate Funktion--- --- --- --- --- --- ----

  # --- --- --- --- --- --- --- --- --- --- --- ---

  # Erstellen der leeren Grundmatrix
  # Erstellen einer Matrix befuellt mit Nullen, die so viel Zeilen und Spalten hat, wie Stellen (liste_stellen) vorhanden sind:
  leere_matrix <- as.data.frame(matrix(data=0, nrow = length(liste_stellen), ncol = length(liste_stellen)))

  # Stellenliste muss nun in Vektor umformatiert werden:
  liste_stellen <- as.vector(liste_stellen)

  # Benennung der Zeilen und Spalten der leeren Matrix
  dimnames(leere_matrix) <- list(liste_stellen, liste_stellen)

  # --- --- --- Befuellen der Matritzen --- --- ---

  # 1. Matrix mit gleich-Werten befuellen:
  matrix_gleich <-befuellen_matrix(leere_matrix, liste_gleich)

  # 2. Matrix mit ueber- und unter-Werten befuellen:
  matrix_ueber_unter <-befuellen_matrix(leere_matrix, liste_ueber)

  # Iteration ueber alle Felder (zunaechst alle Spalten einer Zeile, dann naechste Zeile):
  # Widerspruchsanalyse an Ursprungsdaten ist erforderlich, um kurzkettige Fehler aufzuspueren
  print(paste("Wiederspruchsanalyse laeuft..."))
  for(zeilennr in 1:nrow(matrix_ueber_unter)) # Iteration Zeilen
  {
    for(spaltennr in 1:ncol(matrix_ueber_unter)) # Iteration spalten
    {
      if (matrix_ueber_unter[zeilennr,spaltennr] == 1) # gucken ob Feld besetzt ist in Matrix ueber/unter
      {
        report_of_conflict <- widerspruchsanalyse(zeilennr, spaltennr, TRUE, fehler_loeschen, matrix_gleich, matrix_ueber_unter, liste_ueber_temp, liste_gleich_temp)
        if (report_of_conflict$break_req == TRUE)
        {
        print(report_of_conflict$Conflict)
        print(report_of_conflict$chain_of_conflict)
        return(report_of_conflict)
        }
        else if (report_of_conflict$corr_req == TRUE)
        {
        matrix_ueber_unter <- report_of_conflict$tab_under_above
        matrix_gleich <- report_of_conflict$tab_equal
        }
      }
    }
  }

  # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  # --- --- --- --- --- Vervollstaendigen der Matritzen --- evtl. in separate Funktion--- --- --- --- ----

  # --- --- --- --- --- --- --- --- --- --- --- ---
  durchlaeufe <-1 # fuer die Anzeige, damit man sieht, dass das Programm laeuft
  repeat
  {
    summe_matrix_ueber_unter <- sum(rowSums(matrix_ueber_unter))
    summe_matrix_gleich <- sum(colSums(matrix_gleich))

    # Iteration ueber alle Felder (zunaechst alle Spalten einer Zeile, dann naechste Zeile):
    for(zeilennr in 1:nrow(matrix_ueber_unter)) # Iteration Zeilen
    {
      for(spaltennr in 1:ncol(matrix_ueber_unter)) # Iteration spalten
      {
        if (matrix_gleich[zeilennr,spaltennr] == 1) # gucken ob Feld besetzt ist in Matrix gleich
        {
          matrix_gleich <- gleichsetzen_spalten(matrix_gleich, zeilennr, spaltennr) # Ergaenzen Matrix gleich
          # falls die Diagonale besetzt wird (3=3) wird sie hiermit freigeraeumt, ist nicht notwendig aber effizienter:
          matrix_gleich[zeilennr, zeilennr] <- 0
          # wenn Feld gleich werden Zeilen in Matrix ueber/unter gleichgesetzt:
          matrix_ueber_unter <- gleichsetzen_spalten(matrix_ueber_unter, zeilennr, spaltennr)
          # und auch die Spalten werden gleichgesetzt, daher Transponierung
          matrix_ueber_unter <- t(matrix_ueber_unter)
          matrix_ueber_unter <- gleichsetzen_spalten(matrix_ueber_unter, zeilennr, spaltennr) # Abgleichen Spalten
          matrix_ueber_unter <- t(matrix_ueber_unter) # Matrix wird wieder in Augangslage gekippt
        }

        if (matrix_ueber_unter[zeilennr,spaltennr] == 1) # gucken ob Feld besetzt ist in Matrix ueber/unter
        {
          # Abgleich der Daten innerhalb der "matrix_ueber_unter"
          matrix_ueber_unter <- gleichsetzen_spalten(matrix_ueber_unter, zeilennr, spaltennr) # Gleichsetzen Spalten in der ueber/unter Matrix
          # duerfte Fehler erst nach dem Gleichsetzen der Spalten erkennen, daher vielleicht eine Zeile hoeher? bei Zeiten ausprobieren!
          report_of_conflict <- widerspruchsanalyse(zeilennr, spaltennr, FALSE, fehler_loeschen, matrix_gleich, matrix_ueber_unter, liste_ueber_temp, liste_gleich_temp)
          if (report_of_conflict$break_req == TRUE)
          {
            print(report_of_conflict$Conflict)
            print(report_of_conflict$chain_of_conflict)
            return(report_of_conflict)
          }
          else if (report_of_conflict$corr_req == TRUE)
          {
            matrix_ueber_unter <- report_of_conflict$tab_under_above
            matrix_gleich <- report_of_conflict$tab_equal
          }
        }
      }
    } # Ende Schleife Zeilennr. = 1
    if (summe_matrix_ueber_unter == sum(rowSums(matrix_ueber_unter)) && summe_matrix_gleich == sum(colSums(matrix_gleich))) break()
    print(paste("        Durchlauf:" , durchlaeufe, "abgeschlossen"))
    durchlaeufe <- durchlaeufe +1
  } # Ende Repeat-Schleife

  # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  # --- --- --- --- --- Ergaenzen absoluter Datierungen --- --- --- --- --- --- ----

  # --- --- --- --- --- --- --- --- --- --- --- ---
print("Starting absolute datings...")
  if (missing(absolute_data_list))
  {
    absolute_chronology <- "absolute datings not existing"
  }
  else
  {
    absolute_chronology <- absolute_datierungen(absolute_data_list, liste_stellen, matrix_gleich, matrix_ueber_unter)
     if (absolute_chronology$break_req == TRUE)
    {
      print(absolute_chronology$conflict)
      return(absolute_chronology)
    }
  }
  # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  # --- --- --- --- --- Ausgabe der Daten --- --- --- --- --- --- ----

  # --- --- --- --- --- --- --- --- --- --- --- ---



  stratigraphy <- list(tab_under_above = matrix_ueber_unter, tab_equal = matrix_gleich, absolute_chronology = absolute_chronology$tabelle_dat)

  analysis_functions_objects()

  print("Stratification acomplished")

  return(stratigraphy)


} # End of Function

