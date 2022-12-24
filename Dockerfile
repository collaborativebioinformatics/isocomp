# copied from https://github.com/python-poetry/poetry/issues/1178#issuecomment-998549092
FROM python:3.9-alpine AS builder

WORKDIR /app
ADD pyproject.toml poetry.lock /app/

RUN apk add build-base libffi-dev

RUN pip install --upgrade pip
RUN pip install poetry
RUN poetry config virtualenvs.in-project true
# note --only main means --no-dev. --no-dev is deprecated
RUN poetry install --only main --no-ansi
# alternatively, build from github -- note that you probably need to activate 
# the venv first?
# RUN pip install "git+https://github.com/cmatKhan/pycallingcards.git@raw_processing" --upgrade
# or, better yet, build from pypi

# ---

FROM python:3.9-alpine
WORKDIR /app

COPY --from=builder /app /app
ADD . /app

RUN addgroup -g 1000 app
RUN adduser app -h /app -u 1000 -G 1000 -DH
USER 1000

# change this to match your application
CMD /app/.venv/bin/python -m isocomp