FROM python:3

LABEL maintainer "Hadrien Gourlé <gourlehadrien@gmail.com>"

WORKDIR /usr/src/app

COPY . .
RUN pip install --no-cache-dir .

CMD [ "iss" ]
