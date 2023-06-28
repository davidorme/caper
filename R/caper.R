

# Declare global variables to avoid warnings about usage of variables in
# non-standard evaluation like with() etc.

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "explanatory", "response", "factor", "values", "edge",
        "Nnode", "id", "nTip", "nLin", "tip",
        "N1", "N2", "B", "m", "M", "S.odd", "w", # fusco
        "len", "nSp", "ED", "ED.cor", # pd.calc
    ))
}
