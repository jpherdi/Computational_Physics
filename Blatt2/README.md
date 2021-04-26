# Erklärung

`make run` führt die Datei aus. `make show` führt das Python-Skript aus und zeigt die Plots anschließend im pdf-Viewer der Wahl. Wir nutzen als default `okular` gewählt. Der Viewer kann in der `Makefile` Datei über die Variable `PDF_VIEWER` angepasst werden. `make finish` erzeugt die Pdf-Datei und speichert sie anschließend an die richtige Stelle.

In Latex haben wir einen relativen Pfad eingerichtet, der sich in der `main.tex` ändern lässt. Dieser wird benötigt, damit wir die Latex File mit VSCode automatisiert bauen können, aber auch parallel mit mit der Makefile.