#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <set>
using namespace Rcpp;

// Función auxiliar para encontrar duplicados (equivalente a la lógica de R)
std::vector<int> encontrar_duplicados(const std::vector<int>& ab) {
  std::vector<int> dup;
  std::set<int> vistos;
  std::set<int> duplicados;
  
  for (int elem : ab) {
    if (vistos.count(elem)) {
      duplicados.insert(elem);
    }
    vistos.insert(elem);
  }
  
  for (int elem : ab) {
    if (duplicados.count(elem)) {
      dup.push_back(elem);
    }
  }
  
  return dup;
}

// Función auxiliar para fusionar clusters de manera segura
List fusionar_clusters_seguro(List grupos, std::vector<int>& ab, 
                               const std::vector<int>& B, 
                               const std::vector<int>& dup) {
  if (dup.size() < 2) {
    int A = grupos.size() + (B.size() - ab.size());
    return List::create(Named("grupos") = grupos,
                       Named("ab") = ab,
                       Named("A") = A);
  }
  
  // Buscar en qué clusters está el elemento duplicado
  std::vector<int> ak;
  for (int i = 0; i < grupos.size(); i++) {
    IntegerVector grupo = grupos[i];
    if (std::find(grupo.begin(), grupo.end(), dup[0]) != grupo.end()) {
      ak.push_back(i);
    }
  }
  
  // Validación crítica: verificar que ak tenga suficientes elementos
  if (ak.size() < 2) {
    int A = grupos.size() + (B.size() - ab.size());
    return List::create(Named("grupos") = grupos,
                       Named("ab") = ab,
                       Named("A") = A);
  }
  
  if (dup.size() == 4 && ak.size() >= 2) {
    // Fusionar los clusters que contienen el elemento duplicado
    IntegerVector g1 = grupos[ak[0]];
    IntegerVector g2 = grupos[ak[1]];
    
    // Combinar y hacer único
    std::set<int> unico;
    for (int i = 0; i < g1.size(); i++) unico.insert(g1[i]);
    for (int i = 0; i < g2.size(); i++) unico.insert(g2[i]);
    
    IntegerVector fusionado(unico.size());
    std::copy(unico.begin(), unico.end(), fusionado.begin());
    
    grupos[ak[0]] = fusionado;
    
    // Eliminar el segundo cluster (ajustar índices)
    List grupos_nuevo;
    for (int i = 0; i < grupos.size(); i++) {
      if (i != ak[1]) {
        grupos_nuevo.push_back(grupos[i]);
      }
    }
    grupos = grupos_nuevo;
    
    // Recalcular ab
    ab.clear();
    for (int i = 0; i < grupos.size(); i++) {
      IntegerVector g = grupos[i];
      for (int j = 0; j < g.size(); j++) {
        ab.push_back(g[j]);
      }
    }
    
    // Verificar duplicados nuevamente
    std::vector<int> dup2 = encontrar_duplicados(ab);
    
    // Segunda fusión si aún hay duplicados
    if (dup2.size() >= 2) {
      std::vector<int> ak2;
      for (int i = 0; i < grupos.size(); i++) {
        IntegerVector grupo = grupos[i];
        if (std::find(grupo.begin(), grupo.end(), dup2[0]) != grupo.end()) {
          ak2.push_back(i);
        }
      }
      
      if (ak2.size() >= 2) {
        IntegerVector g1 = grupos[ak2[0]];
        IntegerVector g2 = grupos[ak2[1]];
        
        std::set<int> unico2;
        for (int i = 0; i < g1.size(); i++) unico2.insert(g1[i]);
        for (int i = 0; i < g2.size(); i++) unico2.insert(g2[i]);
        
        IntegerVector fusionado2(unico2.size());
        std::copy(unico2.begin(), unico2.end(), fusionado2.begin());
        
        grupos[ak2[0]] = fusionado2;
        
        List grupos_nuevo2;
        for (int i = 0; i < grupos.size(); i++) {
          if (i != ak2[1]) {
            grupos_nuevo2.push_back(grupos[i]);
          }
        }
        grupos = grupos_nuevo2;
      }
    }
  } else if (ak.size() >= 2) {
    std::vector<int> az = ak;
    
    IntegerVector g1 = grupos[az[0]];
    IntegerVector g2 = grupos[az[az.size() - 1]];
    
    std::set<int> unico;
    for (int i = 0; i < g1.size(); i++) unico.insert(g1[i]);
    for (int i = 0; i < g2.size(); i++) unico.insert(g2[i]);
    
    IntegerVector fusionado(unico.size());
    std::copy(unico.begin(), unico.end(), fusionado.begin());
    
    grupos[az[0]] = fusionado;
    
    List grupos_nuevo;
    for (int i = 0; i < grupos.size(); i++) {
      if (i != az[az.size() - 1]) {
        grupos_nuevo.push_back(grupos[i]);
      }
    }
    grupos = grupos_nuevo;
  }
  
  // Recalcular ab y A
  ab.clear();
  for (int i = 0; i < grupos.size(); i++) {
    IntegerVector g = grupos[i];
    for (int j = 0; j < g.size(); j++) {
      ab.push_back(g[j]);
    }
  }
  
  int A = grupos.size() + (B.size() - ab.size());
  
  return List::create(Named("grupos") = grupos,
                     Named("ab") = ab,
                     Named("A") = A);
}

// [[Rcpp::export]]
List alfa19_cpp(NumericMatrix data, int alpha) {
  int m = data.nrow();
  int n = data.ncol();
  
  // Validación: dataset muy pequeño
  if (m < 3) {
    Rcpp::warning("Dataset muy pequeño. Devolviendo partición trivial.");
    if (alpha == 1) {
      List result;
      IntegerVector todos(m);
      for (int i = 0; i < m; i++) todos[i] = i + 1;
      result.push_back(todos);
      return result;
    }
    if (alpha >= m) {
      List result;
      for (int i = 0; i < m; i++) {
        IntegerVector uno(1);
        uno[0] = i + 1;
        result.push_back(uno);
      }
      return result;
    }
    List result;
    for (int i = 0; i < m; i++) {
      IntegerVector uno(1);
      uno[0] = i + 1;
      result.push_back(uno);
    }
    return result;
  }
  
  // Matriz de coincidencia C
  NumericMatrix C(m, m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      int count = 0;
      for (int k = 0; k < n; k++) {
        if (data(i, k) == data(j, k)) {
          count++;
        }
      }
      C(i, j) = count;
    }
  }
  
  // Matriz de enlace (C sin diagonal)
  NumericMatrix enlace = clone(C);
  for (int i = 0; i < m; i++) {
    enlace(i, i) = 0;
  }
  
  List grupos;
  int A = m;
  std::vector<int> B(m);
  for (int i = 0; i < m; i++) B[i] = i + 1;
  int suma = 1;
  std::vector<int> ab;
  
  // Primer while: Construcción de clusters
  while (A > alpha && suma != 0) {
    std::vector<double> k(m, 0.0);
    for (int j = 0; j < m; j++) {
      double max_val = 0;
      for (int i = 0; i < m; i++) {
        if (enlace(i, j) > max_val) max_val = enlace(i, j);
      }
      k[j] = max_val;
    }
    
    int ini;
    if (A == m) {
      ini = std::max_element(k.begin(), k.end()) - k.begin();
    } else {
      for (int elem : ab) {
        k[elem - 1] = 0;
      }
      double sum_k = 0;
      for (double val : k) sum_k += val;
      if (sum_k == 0) break;
      ini = std::max_element(k.begin(), k.end()) - k.begin();
    }
    
    // Encontrar posición del máximo en la fila ini
    int posicion = 0;
    double max_val = enlace(ini, 0);
    for (int j = 1; j < m; j++) {
      if (enlace(ini, j) > max_val) {
        max_val = enlace(ini, j);
        posicion = j;
      }
    }
    
    double a = enlace(ini, posicion);
    
    // Encontrar todos los índices donde enlace[ini, ] == a
    std::vector<int> di;
    for (int j = 0; j < m; j++) {
      if (enlace(ini, j) == a) {
        di.push_back(j);
      }
    }
    
    enlace(ini, posicion) = 0;
    
    // Encontrar posición del máximo en la columna posicion
    int posicion_columna = 0;
    double max_val_col = enlace(0, posicion);
    for (int i = 1; i < m; i++) {
      if (enlace(i, posicion) > max_val_col) {
        max_val_col = enlace(i, posicion);
        posicion_columna = i;
      }
    }
    
    double b = enlace(posicion_columna, posicion);
    enlace(ini, posicion) = a;
    
    // Procesar según comparación a vs b
    if (a > b || a == b || a < b) {
      if (di.size() > 1) {
        if (m == A) {
          IntegerVector nuevo_grupo(2);
          nuevo_grupo[0] = ini + 1;
          nuevo_grupo[1] = posicion + 1;
          grupos.push_back(nuevo_grupo);
          
          ab.clear();
          for (int i = 0; i < grupos.size(); i++) {
            IntegerVector g = grupos[i];
            for (int j = 0; j < g.size(); j++) {
              ab.push_back(g[j]);
            }
          }
          
          int largo = B.size() - ab.size();
          A = grupos.size() + largo;
          
          for (int i = 0; i < m; i++) {
            enlace(i, posicion) = 0;
            enlace(posicion, i) = 0;
          }
        } else {
          IntegerVector nuevo_grupo(2);
          nuevo_grupo[0] = ini + 1;
          nuevo_grupo[1] = posicion + 1;
          grupos.push_back(nuevo_grupo);
          
          ab.clear();
          for (int i = 0; i < grupos.size(); i++) {
            IntegerVector g = grupos[i];
            for (int j = 0; j < g.size(); j++) {
              ab.push_back(g[j]);
            }
          }
          
          std::vector<int> dup = encontrar_duplicados(ab);
          List resultado = fusionar_clusters_seguro(grupos, ab, B, dup);
          grupos = resultado["grupos"];
          ab = as<std::vector<int>>(resultado["ab"]);
          A = as<int>(resultado["A"]);
          
          for (int i = 0; i < m; i++) {
            enlace(i, posicion) = 0;
            enlace(posicion, i) = 0;
          }
        }
        
        // Procesar di2
        std::vector<int> di2;
        for (int j = 0; j < m; j++) {
          if (enlace(ini, j) == a) {
            di2.push_back(j);
          }
        }
        
        for (size_t jj = 0; jj < di2.size(); jj++) {
          if (A > alpha) {
            // Encontrar posición del máximo en la fila ini
            int pos2 = 0;
            double max_val2 = enlace(ini, 0);
            for (int j = 1; j < m; j++) {
              if (enlace(ini, j) > max_val2) {
                max_val2 = enlace(ini, j);
                pos2 = j;
              }
            }
            
            double a2 = enlace(ini, pos2);
            enlace(ini, pos2) = 0;
            
            // Encontrar posición del máximo en la columna pos2
            int pos_col2 = 0;
            double max_val_col2 = enlace(0, pos2);
            for (int i = 1; i < m; i++) {
              if (enlace(i, pos2) > max_val_col2) {
                max_val_col2 = enlace(i, pos2);
                pos_col2 = i;
              }
            }
            
            double b2 = enlace(pos_col2, pos2);
            enlace(ini, pos2) = a2;
            
            if (a2 >= b2) {
              for (int i = 0; i < m; i++) {
                enlace(i, pos2) = 0;
                enlace(pos2, i) = 0;
              }
              
              IntegerVector nuevo_grupo2(2);
              nuevo_grupo2[0] = ini + 1;
              nuevo_grupo2[1] = pos2 + 1;
              grupos.push_back(nuevo_grupo2);
              
              ab.clear();
              for (int i = 0; i < grupos.size(); i++) {
                IntegerVector g = grupos[i];
                for (int j = 0; j < g.size(); j++) {
                  ab.push_back(g[j]);
                }
              }
              
              std::vector<int> dup = encontrar_duplicados(ab);
              List resultado = fusionar_clusters_seguro(grupos, ab, B, dup);
              grupos = resultado["grupos"];
              ab = as<std::vector<int>>(resultado["ab"]);
              A = as<int>(resultado["A"]);
            }
          }
        }
      } else {
        if (m == A) {
          IntegerVector nuevo_grupo(2);
          nuevo_grupo[0] = ini + 1;
          nuevo_grupo[1] = posicion + 1;
          grupos.push_back(nuevo_grupo);
          
          ab.clear();
          for (int i = 0; i < grupos.size(); i++) {
            IntegerVector g = grupos[i];
            for (int j = 0; j < g.size(); j++) {
              ab.push_back(g[j]);
            }
          }
          
          int largo = B.size() - ab.size();
          A = grupos.size() + largo;
          
          for (int i = 0; i < m; i++) {
            enlace(i, posicion) = 0;
            enlace(posicion, i) = 0;
          }
        } else {
          IntegerVector nuevo_grupo(2);
          nuevo_grupo[0] = ini + 1;
          nuevo_grupo[1] = posicion + 1;
          grupos.push_back(nuevo_grupo);
          
          for (int i = 0; i < m; i++) {
            enlace(i, posicion) = 0;
            enlace(posicion, i) = 0;
          }
          
          ab.clear();
          for (int i = 0; i < grupos.size(); i++) {
            IntegerVector g = grupos[i];
            for (int j = 0; j < g.size(); j++) {
              ab.push_back(g[j]);
            }
          }
          
          std::vector<int> dup = encontrar_duplicados(ab);
          List resultado = fusionar_clusters_seguro(grupos, ab, B, dup);
          grupos = resultado["grupos"];
          ab = as<std::vector<int>>(resultado["ab"]);
          A = as<int>(resultado["A"]);
        }
      }
    }
    
    // Actualizar k y verificar terminación
    std::fill(k.begin(), k.end(), 0.0);
    for (int j = 0; j < m; j++) {
      double max_val = 0;
      for (int i = 0; i < m; i++) {
        if (enlace(i, j) > max_val) max_val = enlace(i, j);
      }
      k[j] = max_val;
    }
    
    ab.clear();
    for (int i = 0; i < grupos.size(); i++) {
      IntegerVector g = grupos[i];
      for (int j = 0; j < g.size(); j++) {
        ab.push_back(g[j]);
      }
    }
    
    for (int elem : ab) {
      k[elem - 1] = 0;
    }
    
    int largo = B.size() - ab.size();
    
    double sum_k = 0;
    for (double val : k) sum_k += val;
    
    if (sum_k == 0 || largo == 0) {
      suma = 0;
    } else {
      suma = sum_k;
    }
    
    int k2 = 0;
    for (double val : k) {
      if (val != 0) k2++;
    }
    
    if (k2 == 1) {
      suma = 0;
    }
  }
  
  // Segundo while: Optimización de clusters
  while (A > alpha && suma == 0) {
    int ab3 = 0;
    for (int i = 0; i < grupos.size(); i++) {
      IntegerVector g = grupos[i];
      ab3 += g.size();
    }
    
    if (ab3 == m) {
      int num_grupos = grupos.size();
      int num_comb = (num_grupos * (num_grupos - 1)) / 2;
      
      std::vector<std::pair<int, int>> combinaciones;
      for (int i = 0; i < num_grupos; i++) {
        for (int j = i + 1; j < num_grupos; j++) {
          combinaciones.push_back(std::make_pair(i, j));
        }
      }
      
      std::vector<double> Calidad_c(num_comb, 0.0);
      
      for (int kp = 0; kp < num_comb; kp++) {
        List grupos2 = clone(grupos);
        
        IntegerVector g1 = grupos2[combinaciones[kp].first];
        IntegerVector g2 = grupos2[combinaciones[kp].second];
        
        IntegerVector fusionado(g1.size() + g2.size());
        for (int i = 0; i < g1.size(); i++) fusionado[i] = g1[i];
        for (int i = 0; i < g2.size(); i++) fusionado[g1.size() + i] = g2[i];
        
        grupos2[combinaciones[kp].first] = fusionado;
        
        // Eliminar el segundo grupo
        List grupos2_final;
        for (int i = 0; i < grupos2.size(); i++) {
          if (i != combinaciones[kp].second) {
            grupos2_final.push_back(grupos2[i]);
          }
        }
        
        // Crear class_C
        std::vector<int> class_C(m, 0);
        for (int g = 0; g < grupos2_final.size(); g++) {
          IntegerVector grupo = grupos2_final[g];
          for (int i = 0; i < grupo.size(); i++) {
            class_C[grupo[i] - 1] = g + 1;
          }
        }
        
        // Calcular S y QI
        double termino1 = 0.0;
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            double s_ij = (class_C[i] == class_C[j]) ? 1.0 : 0.0;
            termino1 += (3.0 * C(i, j) - 2.0 * n) * s_ij;
          }
        }
        
        double termino2 = 0.0;
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            termino2 += 2.0 * (n - C(i, j));
          }
        }
        
        double numerador = termino1 + termino2;
        double denominador = (double)m * m * n;
        Calidad_c[kp] = numerador / denominador;
      }
      
      int max_calidad = std::max_element(Calidad_c.begin(), Calidad_c.end()) - Calidad_c.begin();
      int idx1 = combinaciones[max_calidad].first;
      int idx2 = combinaciones[max_calidad].second;
      
      IntegerVector g1 = grupos[idx1];
      IntegerVector g2 = grupos[idx2];
      
      IntegerVector fusionado(g1.size() + g2.size());
      for (int i = 0; i < g1.size(); i++) fusionado[i] = g1[i];
      for (int i = 0; i < g2.size(); i++) fusionado[g1.size() + i] = g2[i];
      
      grupos[idx1] = fusionado;
      
      List grupos_nuevo;
      for (int i = 0; i < grupos.size(); i++) {
        if (i != idx2) {
          grupos_nuevo.push_back(grupos[i]);
        }
      }
      grupos = grupos_nuevo;
      
      ab.clear();
      for (int i = 0; i < grupos.size(); i++) {
        IntegerVector g = grupos[i];
        for (int j = 0; j < g.size(); j++) {
          ab.push_back(g[j]);
        }
      }
      
      int largo = B.size() - ab.size();
      A = grupos.size() + largo;
    } else {
      ab.clear();
      for (int i = 0; i < grupos.size(); i++) {
        IntegerVector g = grupos[i];
        for (int j = 0; j < g.size(); j++) {
          ab.push_back(g[j]);
        }
      }
      
      std::vector<int> largo2;
      for (int elem : B) {
        if (std::find(ab.begin(), ab.end(), elem) == ab.end()) {
          largo2.push_back(elem);
        }
      }
      
      IntegerVector nuevo_grupo(largo2.size());
      for (size_t i = 0; i < largo2.size(); i++) {
        nuevo_grupo[i] = largo2[i];
      }
      grupos.push_back(nuevo_grupo);
      
      int num_grupos = grupos.size();
      int num_comb = (num_grupos * (num_grupos - 1)) / 2;
      
      std::vector<std::pair<int, int>> combinaciones;
      for (int i = 0; i < num_grupos; i++) {
        for (int j = i + 1; j < num_grupos; j++) {
          combinaciones.push_back(std::make_pair(i, j));
        }
      }
      
      std::vector<double> Calidad_c(num_comb, 0.0);
      
      for (int kp = 0; kp < num_comb; kp++) {
        List grupos2 = clone(grupos);
        
        IntegerVector g1 = grupos2[combinaciones[kp].first];
        IntegerVector g2 = grupos2[combinaciones[kp].second];
        
        IntegerVector fusionado(g1.size() + g2.size());
        for (int i = 0; i < g1.size(); i++) fusionado[i] = g1[i];
        for (int i = 0; i < g2.size(); i++) fusionado[g1.size() + i] = g2[i];
        
        grupos2[combinaciones[kp].first] = fusionado;
        
        List grupos2_final;
        for (int i = 0; i < grupos2.size(); i++) {
          if (i != combinaciones[kp].second) {
            grupos2_final.push_back(grupos2[i]);
          }
        }
        
        std::vector<int> class_C(m, 0);
        for (int g = 0; g < grupos2_final.size(); g++) {
          IntegerVector grupo = grupos2_final[g];
          for (int i = 0; i < grupo.size(); i++) {
            class_C[grupo[i] - 1] = g + 1;
          }
        }
        
        double termino1 = 0.0;
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            double s_ij = (class_C[i] == class_C[j]) ? 1.0 : 0.0;
            termino1 += (3.0 * C(i, j) - 2.0 * n) * s_ij;
          }
        }
        
        double termino2 = 0.0;
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            termino2 += 2.0 * (n - C(i, j));
          }
        }
        
        double numerador = termino1 + termino2;
        double denominador = (double)m * m * n;
        Calidad_c[kp] = numerador / denominador;
      }
      
      int max_calidad = std::max_element(Calidad_c.begin(), Calidad_c.end()) - Calidad_c.begin();
      int idx1 = combinaciones[max_calidad].first;
      int idx2 = combinaciones[max_calidad].second;
      
      IntegerVector g1 = grupos[idx1];
      IntegerVector g2 = grupos[idx2];
      
      IntegerVector fusionado(g1.size() + g2.size());
      for (int i = 0; i < g1.size(); i++) fusionado[i] = g1[i];
      for (int i = 0; i < g2.size(); i++) fusionado[g1.size() + i] = g2[i];
      
      grupos[idx1] = fusionado;
      
      List grupos_nuevo;
      for (int i = 0; i < grupos.size(); i++) {
        if (i != idx2) {
          grupos_nuevo.push_back(grupos[i]);
        }
      }
      grupos = grupos_nuevo;
      
      ab.clear();
      for (int i = 0; i < grupos.size(); i++) {
        IntegerVector g = grupos[i];
        for (int j = 0; j < g.size(); j++) {
          ab.push_back(g[j]);
        }
      }
      
      int largo = B.size() - ab.size();
      A = grupos.size() + largo;
    }
  }
  
  // Asignación final
  List z = clone(grupos);
  
  ab.clear();
  for (int i = 0; i < grupos.size(); i++) {
    IntegerVector g = grupos[i];
    for (int j = 0; j < g.size(); j++) {
      ab.push_back(g[j]);
    }
  }
  
  std::vector<int> ac;
  for (int elem : B) {
    if (std::find(ab.begin(), ab.end(), elem) == ab.end()) {
      ac.push_back(elem);
    }
  }
  
  if ((int)grupos.size() == alpha) {
    z = grupos;
  } else {
    for (size_t i = 0; i < ac.size(); i++) {
      IntegerVector uno(1);
      uno[0] = ac[i];
      z.push_back(uno);
    }
  }
  
  return z;
}