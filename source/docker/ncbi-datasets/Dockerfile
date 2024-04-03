FROM alpine:3.18.3

# Install bash and curl
RUN apk add --no-cache bash curl

# Install ncbi-datasets
RUN curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets' && \
    curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat' && \
    chmod +x datasets dataformat && \
    mv datasets dataformat /usr/local/bin/

# set entrypoint
CMD [ "bash" ]