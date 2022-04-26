PNG.vec <- Sys.glob("figure-aum-convexity-interactive-screenshots/*.PNG")
for(PNG in PNG.vec){
  img <- magick::image_read(PNG)
  crop <- magick::image_crop(img, "1455x950+17+93")
  out <- sub("screenshots", "cropped", PNG)
  dir.create(dirname(out), showWarnings=FALSE, recursive=TRUE)
  magick::image_write(crop, out)
cat(sprintf(r"{
\begin{frame}
  \frametitle{Demonstration of AUC/AUM computation}
  {\scriptsize\url{https://bl.ocks.org/tdhock/raw/545d76ea8c0678785896e7dbe5ff5510/}}

  \includegraphics[width=\textwidth]{%s}
\end{frame}
}", out))
}
