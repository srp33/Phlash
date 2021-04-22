import 'bootstrap/dist/css/bootstrap.css';
import 'bootstrap-vue/dist/bootstrap-vue.css';
import BootstrapVue from 'bootstrap-vue';
import Vue from 'vue';
import App from './App.vue';
import router from './router';
import GoogleLogin from 'vue-google-login';
import { LoaderPlugin } from 'vue-google-login';

Vue.use(BootstrapVue);
Vue.use(LoaderPlugin, {
  client_id: process.env.VUE_APP_API_KEY
});

Vue.GoogleAuth.then(auth2 => {
  if (!auth2.isSignedIn.get()) {
    router.push('/');
  }
})
Vue.config.productionTip = false;

new Vue({
  router,
  GoogleLogin,
  LoaderPlugin,
  render: (h) => h(App),
}).$mount('#app');
