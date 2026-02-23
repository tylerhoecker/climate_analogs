# Climate Analogs

**Analyses to identify climate analogs and implement climate-analog impact models**

This repository contains dual implementations (R and Julia) for identifying climate analogs using Mahalanobis distance and applying them to predict vegetation shifts under climate change.


---

## Table of Contents

* [Overview](#overview)
* [Project Workflow](#project-workflow)
* [File Structure](#file-structure)
* [Key Variables & Concepts](#key-variables--concepts)
* [Core Functions](#core-functions)
* [Dependencies](#dependencies)
* [Usage Examples](#usage-examples)


---

## Overview

**Reverse climate analogs** are geographic locations that currently experience climate conditions similar to what a focal location will experience under future climate change. Forward analogs are the opposite. This project:


1. **Calculates climate dissimilarity** from a focal location and its future climate to the contemporary climate mean of surrounding analog candidates using Mahalanobis distance (MD) transformed to sigma (σ) values
2. **Identifies the best analogs** from a historical climate pool
3. **Predicts future impacts** (e.g., vegetation shifts) by extracting vegetation at analog locations

   ## Project Workflow

 ![Flowchart of the climate analogs process performed](./flowchart.svg)

## File Structure

 ![File Structure](file_structure.svg)