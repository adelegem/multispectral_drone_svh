# Pin to the R version captured in renv.lock so the lockfile is satisfiable.
# rocker/geospatial ships GDAL, PROJ, GEOS, and UDUNITS pre-installed —
# saves us the system-dep dance that GitHub Actions has to do.
FROM rocker/geospatial:4.5.3

LABEL org.opencontainers.image.source="https://github.com/adelegem/multispectral_drone_svh"
LABEL org.opencontainers.image.description="Reproducible analysis environment for Gemmell et al., applying the spectral variability hypothesis to arid shrublands"
LABEL org.opencontainers.image.licenses="MIT"

# Optional GitHub PAT for renv to pull saltbush (avoids rate limiting).
# Pass at build time with `--build-arg GITHUB_PAT=$(gh auth token)`.
ARG GITHUB_PAT
ENV GITHUB_PAT=${GITHUB_PAT}

WORKDIR /workspace

# System libraries not bundled with rocker/geospatial: glpk + gmp for igraph
# (a ggraph dependency). Everything else (sf, terra, units, ragg, textshaping,
# xml2, curl, openssl, ...) is supplied by the base image.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libglpk-dev libgmp-dev \
    && rm -rf /var/lib/apt/lists/*

# Restore the locked R environment first so the heavy package-install layer
# is cached across changes to anything else in the project.
COPY .Rprofile renv.lock ./
COPY renv/activate.R renv/settings.json renv/

RUN R -e "renv::restore(prompt = FALSE)"

# Now copy the rest of the project (everything not excluded by .dockerignore).
COPY . .

# Default command runs the full pipeline. Override at runtime with e.g.
# `docker run --rm -it ... R` to drop into an interactive session.
CMD ["Rscript", "-e", "targets::tar_make()"]
