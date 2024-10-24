#!/bin/sh

workdir=${HOME}/rstudio-tmp

mkdir -p -m 700 ${workdir}/run ${workdir}/tmp ${workdir}/var/lib/rstudio-server
cat > ${workdir}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

# Set R_LIBS_USER to a specific path to avoid conflicts with
# personal libraries from any R installation in the host environment
cat > ${workdir}/rsession.sh <<END
#!/bin/sh
export R_LIBS_USER=${HOME}/R/my_rstudio/4.4.0
exec rsession "\${@}"
END

chmod +x ${workdir}/rsession.sh

# Do not suspend idle sessions.
# Alternative to setting session-timeout-minutes=0 in /etc/rstudio/rsession.conf
# https://github.com/rstudio/rstudio/blob/v1.4.1106/src/cpp/server/ServerSessionManager.cpp#L126
export SINGULARITYENV_RSTUDIO_SESSION_TIMEOUT=0

export SINGULARITYENV_USER=$(id -un)

# Get unused socket per https://unix.stackexchange.com/a/132524
# Tiny race condition between the python & singularity commands
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
# Get node IP address.
readonly ADD=$(nslookup `hostname` | grep -i address | awk -F" " '{print $2}' | awk -F# '{print $1}' | tail -n 1)

cat 1>&2 <<END
"Running RStudio at $ADD:$PORT"
END

singularity exec --cleanenv -c -W ${workdir} --bind ${HOME}:/home/$USER,/research:/research,/hpcf:/hpcf,${workdir}/run:/run,${workdir}/tmp:/tmp,${workdir}/database.conf:/etc/rstudio/database.conf,${workdir}/rsession.sh:/etc/rstudio/rsession.sh,${workdir}/var/lib/rstudio-server:/var/lib/rstudio-server rstudio_4.4.0.sif \
    rserver --www-port ${PORT} \
            --rsession-path=/etc/rstudio/rsession.sh \
            --secure-cookie-key-file ${workdir}/tmp/rstudio-server/secure-cookie-key \
            --auth-stay-signed-in-days=30 \
            --auth-none=1 \
            --server-user $USER \
            --auth-timeout-minutes=0
