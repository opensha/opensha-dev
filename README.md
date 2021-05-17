# opensha-dev

Development sandbox for OpenSHA project

This project depends on the main [opensha repository](https://github.com/opensha/opensha), and is intended for exporatory and test code. Once you have checkout of the `opensha` project, follow the following steps to set up this project:

### Cloning in a terminal

To clone this project in a terminal (simplest and quickest for most users):

```bash
cd git # or whatever directory you previously checked out the 'opensha' repository into
git clone https://github.com/opensha/opensha-dev.git
```

### Building in a terminal with Gradle

OpenSHA uses Gradle to handle the build process from dependency management to compilation. You should compile OpenSHA from your lowest level project, *opensha-dev* for the example above.

```bash
cd opensha-dev
./gradlew assemble
```

This will build all source files in opensha-dev and parent projects. It will also build a jar file for each project, not including any dependencies. You can build a "fat jar" which includes dependent libraries as follows:

```bash
cd opensha-dev # or whichever project you are interested in
./gradlew fatJar
```

### Developing & building OpenSHA with Eclipse

Most active OpenSHA development is done through [Eclipse](https://eclipse.org). You will need the Eclipse IDE for Java Developers.

>**NOTE:** The following instructions assume that you have already cloned the OpenSHA projects on a terminal, though you can clone them through Eclipse. If you chose to go this route, be sure to leave the "Import all existing Eclipse projects after clone finishes" check-box **UNSELECTED**, as this feature will cause issues with gradle.

For each project, you will need to do the following:
* `File > Import`  
* Select `Gradle > Existing Gradle Project` and hit `Next`  
* Browse to the location of `opensha` under `Project root directory`  
* Hit `Finish`  
* Repeat for this (and any other) sub-projects. **IMPORTANT: projects must be imported in order, dependent projects first. That means opensha, then opensha-dev**  

You can either use Eclipse's built in Git tools, or the Git command line client to pull/push changes. If any of the `.gradle` files are modified, right click on the project within eclipse and select `Gradle >  Refresh Gradle Project`
