FROM node:lts-alpine as build-stage

WORKDIR /app

COPY package*.json ./
RUN npm install
#RUN npm audit fix
COPY ./ .
COPY .env.production .
RUN npm run build

FROM nginx as production-stage
RUN mkdir /app
COPY --from=build-stage /app/dist /app/phlash
RUN mv /app/phlash/index.html /app/
COPY nginx.conf /etc/nginx/nginx.conf

EXPOSE 80
