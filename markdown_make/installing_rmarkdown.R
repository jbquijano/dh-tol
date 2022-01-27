# This R script includes the instructions for
# installing rmarkdown, bookdown and other markdown
# extension packages

# R markdown===============================================

# Install rmarkdown
install.packages('rmarkdown')
# I like using the pacman package in installing 
# and updating stuff because IT INSTALLS AND UPDATES
# (installs if not yet installed and update if already
# installed but there are package updates)
pacman::p_load('rmarkdown')

# Do you want PDF outputs? Yep! Then install LaTeX.
# If not yet installed, TinyTex can be installed using R
install.packages('tinytex')
tinytex::install_tinytex()  # install TinyTeX

# bookdown================================================
# Rmarkdown is fine, but if you want a book-like appearance
# to your markdown document, you can use bookdown.
# This looks preetier and more professional (I like this
# for its better table of contents and paged html)

# Install bookdown
# static version
install.packages("bookdown")
# or development version on GitHub
remotes::install_github('rstudio/bookdown')
# you can also install through pacman
pacman::p_load(bookdown)

# klippy===================================================

# This package adds a copy to clipboard buttons in
# R markdown html docs.

# Install klippy
remotes::install_github("rlesur/klippy")
