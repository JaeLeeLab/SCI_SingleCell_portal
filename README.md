
## Shiny Single Cell Browser

## Deploy the Single Cell Browser
  - The App can be easily deployed on a web server using [shinyapps.io](https://www.shinyapps.io), which supports both free and paid servers. Docker is an alternative approach but takes longer to set up.
  - To set up a [shinyapps.io](https://www.shinyapps.io) account and learn how to deploy a Shiny app, follow this [tutorial](https://shiny.rstudio.com/articles/shinyapps.html).
  - After setting up the account, deploy the app by [`rsconnect::deployApp()`](https://www.rdocumentation.org/packages/rsconnect/versions/0.8.18/topics/deployApp).

If you encounter the following error: `Error parsing manifest: Unable to determine package source for Bioconductor package Biobase: Repository must be specified`, run this before deployApp: `options(repos = BiocManager::repositories()`

  
## Updates

See [updates.md](UPDATES.md)



