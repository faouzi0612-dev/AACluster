#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix calcular_matriz_C_cpp(NumericMatrix data) {
  int m = data.nrow();
  int n = data.ncol();
  NumericMatrix C(m, m);
  
  // Calcular matriz de similitud C (coincidencia exacta)
  for(int i = 0; i < m; i++) {
    for(int j = i; j < m; j++) {
      int count = 0;
      for(int k = 0; k < n; k++) {
        if(data(i, k) == data(j, k)) {
          count++;
        }
      }
      C(i, j) = count;
      C(j, i) = count; // Simétrica
    }
  }
  
  return C;
}

// [[Rcpp::export]]
NumericMatrix calcular_matriz_Y_c_cpp(IntegerVector class_C, int m) {
  NumericMatrix Y_c(m, m);
  
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < m; j++) {
      Y_c(i, j) = (class_C[i] == class_C[j]) ? 1.0 : 0.0;
    }
  }
  
  return Y_c;
}

// [[Rcpp::export]]
NumericVector calcular_calidad_c_cpp(NumericMatrix C, List combinaciones, 
                                     List grupos, int m, int n) {
  int num_comb = combinaciones.size();
  NumericVector Calidad_c(num_comb);
  
  // Pre-calcular denominador (constante, no cambia)
   double denominador = (2.0 - 1.0 / m) * m * m * n;
    // double denominador =   n*(m^2 +2*m-2);
  
  // Pre-calcular termino2 (constante, no depende de la partición)
  double termino2 = 0.0;
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < m; j++) {
      termino2 += 2.0 * (n - C(i, j));
    }
  }
  
  for(int kp = 0; kp < num_comb; kp++) {
    IntegerVector comb = combinaciones[kp];
    List grupos2 = clone(grupos);
    IntegerVector idx1 = comb[0] - 1;
    IntegerVector idx2 = comb[1] - 1;
    
    List g1 = grupos2[idx1[0]];
    List g2 = grupos2[idx2[0]];
    
    IntegerVector fusionado(g1.size() + g2.size());
    for(int i = 0; i < g1.size(); i++) fusionado[i] = g1[i];
    for(int i = 0; i < g2.size(); i++) fusionado[g1.size() + i] = g2[i];
    
    grupos2[idx1[0]] = fusionado;
    grupos2[idx2[0]] = R_NilValue;
    
    List grupos2_final;
    for(int i = 0; i < grupos2.size(); i++) {
      if(!Rf_isNull(grupos2[i])) {
        grupos2_final.push_back(grupos2[i]);
      }
    }
    
    IntegerVector class_C(m);
    for(int g = 0; g < grupos2_final.size(); g++) {
      IntegerVector grupo = grupos2_final[g];
      for(int i = 0; i < grupo.size(); i++) {
        class_C[grupo[i] - 1] = g + 1;
      }
    }
    
    // Calcular termino1: sum((3*C - 2*n) * S)
    double termino1 = 0.0;
    for(int i = 0; i < m; i++) {
      for(int j = 0; j < m; j++) {
        double s_ij = (class_C[i] == class_C[j]) ? 1.0 : 0.0;
        termino1 += (3.0 * C(i, j) - 2.0 * n) * s_ij;
      }
    }
    
    // Fórmula correcta del QI
    double numerador = termino1 + termino2;
    Calidad_c[kp] = numerador / denominador;
  }
  
  return Calidad_c;
}
// [[Rcpp::export]]
double calcular_QI_particion_cpp(NumericMatrix C, List grupos, int m, int n) {
  // Crear matriz de membresía S
  NumericMatrix S(m, m);
  for(int k = 0; k < grupos.size(); k++) {
    IntegerVector idx = grupos[k];
    for(int i = 0; i < idx.size(); i++) {
      for(int j = 0; j < idx.size(); j++) {
        S(idx[i] - 1, idx[j] - 1) = 1.0;
      }
    }
  }
  
  // Calcular términos
  double termino1 = 0.0;
  double termino2 = 0.0;
  
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < m; j++) {
      termino1 += (3.0 * C(i, j) - 2.0 * n) * S(i, j);
      termino2 += 2.0 * (n - C(i, j));
    }
  }
  
  double numerador = termino1 + termino2;
  double denominador = (2.0 - 1.0 / m) * m * m * n;
  
  return numerador / denominador;
}