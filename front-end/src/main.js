import 'bootstrap/dist/css/bootstrap.css';
import 'bootstrap-vue/dist/bootstrap-vue.css';
import BootstrapVue from 'bootstrap-vue';
import Vue from 'vue';
import App from './App.vue';
import router from './router';
import GoogleSignInButton from 'vue-google-signin-button-directive';
import GSignInButton from 'vue-google-signin-button';

Vue.use(BootstrapVue);
Vue.use(GSignInButton);
Vue.config.productionTip = false;

new Vue({
  router,
  GoogleSignInButton,
  GSignInButton,
  render: (h) => h(App),
}).$mount('#app');
