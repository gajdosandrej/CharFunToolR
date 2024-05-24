plotReIm3 <- function(f, x, title = "", xmin, xmax, dx = 0.1, xlab = "t", ylab = "f(t)",
                      px = 0, py = 0, labelx = TRUE, labely = TRUE) {
        if(missing(xmin)) {
                xmin <- vector()
        }
        if(missing(xmax)) {
                xmax <- vector()
        }


        f_real <- f
        f_imag <- f



        for(i in 1:length(f)){

                for (h in 1:length(x)) {


                f_real[[i]][h] <- Re(f[[i]][h])
                f_imag[[i]][h] <- Im(f[[i]][h])
        }
        }
        plot(NULL, main = title, xlab = xlab, ylab = ylab, xlim = c(min(x), max(x)),
             ylim = c(min(c(unlist(f_real), unlist(f_imag))), max(c(unlist(f_real), unlist(f_imag)))))

      #  if (labelx) {
       #         for (xi in c(px)) {
        #                lines(c(xi, xi),
         #                     c(min(c(unlist(f_real), unlist(f_imag))) - 1, max(c(unlist(f_real), unlist(f_imag))) + 1),
          #                    lty = 1,
            #                  lwd = 2,
           #                   col = "gray"
                   #     )
             #   }
        #}
        #if (labely) {
         #       for (yi in c(py)) {
          #              lines(
           #                     c(min(x), max(x)),
            #                    c(yi, yi),
             #                   lty = 1,
              #                  lwd = 2,
               #                 col = "gray"
                #        )
                #}
        #}

        lines_col <- randomcoloR::distinctColorPalette(2 * length(f))
        for (j in 1 : length(f)) {



                lines(x, f_real[[j]], lty = 1, lwd = 2, col = lines_col[2 * j - 1])
                lines(x, f_imag[[j]], lty = 1, lwd = 2, col = lines_col[2 * j])
        }

        legend("bottomright",
               legend = c("Re(f))",rep("", length(f) - 1),"Im(f)",rep("", length(f) - 1)),
               fill = lines_col,
               bg = "white", ncol = 2, cex = 0.75)
        grid()

}
