import math
from scipy import stats

def calculate_independent_t_test_power(effect_size, n1, n2, alpha=0.05):
    """
    G*Power mantığıyla bağımsız örneklemler t-testi için A priori güç (power) analizi yapar.
    """
    # Serbestlik derecesi
    df = n1 + n2 - 2
    
    # Kritik t değeri (çift yönlü test için)
    t_crit = stats.t.ppf(1 - alpha / 2, df)
    
    # Merkezi olmayan t dağılımı için non-centrality parameter (delta)
    ncp = effect_size * math.sqrt((n1 * n2) / (n1 + n2))
    
    # İstatistiksel güç (1 - beta) hesaplaması
    # Güç = P(T > t_crit | ncp) + P(T < -t_crit | ncp)
    power = 1 - (stats.nct.cdf(t_crit, df, ncp) - stats.nct.cdf(-t_crit, df, ncp))
    
    return t_crit, ncp, power

def calculate_required_sample_size_cohen_d(effect_size, alpha=0.05, target_power=0.80):
    """
    Hedeflenen istatistiksel güce (örn: 0.80) ulaşmak için gereken örneklem boyutunu bulur.
    """
    n = 2  # Başlangıç örneklem boyutu
    while True:
        df = 2 * n - 2
        t_crit = stats.t.ppf(1 - alpha / 2, df)
        ncp = effect_size * math.sqrt((n * n) / (2 * n))
        power = 1 - (stats.nct.cdf(t_crit, df, ncp) - stats.nct.cdf(-t_crit, df, ncp))
        
        if power >= target_power:
            return n, power
        n += 1

if __name__ == "__main__":
    print("--- G*Power Simülasyon Motoru ---")
    
    # Örnek Senaryo: Cohen's d = 0.5 (Orta etki büyüklüğü), her grup için n=30
    es = 0.5
    group1_n = 30
    group2_n = 30
    
    t_val, non_centrality, calculated_power = calculate_independent_t_test_power(es, group1_n, group2_n)
    print(f"Etki Büyüklüğü (d): {es}")
    print(f"Grup Örneklemleri: n1={group1_n}, n2={group2_n}")
    print(f"Hesaplanan İstatistiksel Güç (1-beta): {calculated_power:.4f}")
    
    # Örnek Senaryo: %80 güç için gereken örneklem boyutu
    req_n, achieved_power = calculate_required_sample_size_cohen_d(effect_size=0.5, target_power=0.80)
    print("\n--- A priori Örneklem Boyutu Analizi ---")
    print(f"Hedef Güç: 0.80 | Etki Büyüklüğü: 0.5")
    print(f"Gerekli Tek Grup Örneklem Boyutu (n): {req_n} (Toplam: {req_n * 2})")
    print(f"Elde Edilen Güç: {achieved_power:.4f}")