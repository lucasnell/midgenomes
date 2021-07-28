
if (interactive()) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 8, height = 6,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}
