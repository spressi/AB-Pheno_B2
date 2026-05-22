allCols = "Time	Type	Trial	L Raw X [px]	L Raw Y [px]	R Raw X [px]	R Raw Y [px]	L Dia X [px]	L Dia Y [px]	L Mapped Diameter [mm]	R Dia X [px]	R Dia Y [px]	R Mapped Diameter [mm]	L CR1 X [px]	L CR1 Y [px]	L CR2 X [px]	L CR2 Y [px]	R CR1 X [px]	R CR1 Y [px]	R CR2 X [px]	R CR2 Y [px]	L POR X [px]	L POR Y [px]	R POR X [px]	R POR Y [px]	Timing	L Validity	R Validity	Pupil Confidence	L Plane	R Plane	L EPOS X	L EPOS Y	L EPOS Z	R EPOS X	R EPOS Y	R EPOS Z	L GVEC X	L GVEC Y	L GVEC Z	R GVEC X	R GVEC Y	R GVEC Z	Trigger	Frame	Aux1" %>% 
  strsplit("\t") %>% unlist()

targetCols = "Time	Type	Trial	L Dia X [px]	L Dia Y [px]	L Mapped Diameter [mm]	R Dia X [px]	R Dia Y [px]	R Mapped Diameter [mm]	L POR X [px]	L POR Y [px]	R POR X [px]	R POR Y [px]" %>% 
  strsplit("\t") %>% unlist()

targetIndeces = which(allCols %in% targetCols)
