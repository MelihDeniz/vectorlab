"""PKU vektörünün insanda beklenen etkisini kaba bir modelle tahmin eder.

Bu bir TAHMİN ARACIDIR, bir kehanet değil. Modelin amacı "bizim vektör
şu kadar düşürür" demek değil; **ne kadar aktivite gerekirse ne olur**
sorusunu görünür kılmaktır.

Model tek bölmeli bir kütle dengesidir:

    dC/dt = R_ogun(t) + k_salim - k_sentez * C - D(t, C)

    C          plazma fenilalanin derişimi (µM)
    R_ogun(t)  öğünden emilen Phe (µM/sa), günde üç öğün
    k_salim    doku proteini yıkımından sürekli Phe salımı (µM/sa)
    k_sentez   protein sentezine geri giden Phe (1/sa)
    D(t, C)    bakterinin bağırsakta tükettiği Phe (µM/sa),
               yalnızca öğün sonrası geçiş penceresinde, doyuma giden
               (Michaelis-Menten) biçimde

PKU'da fenilalanin hidroksilaz çalışmadığı için modelde o yol yoktur:
tedavisiz durumda tek çıkış protein sentezidir, bu yüzden Phe birikir.

KALİBRASYON. Bakteri terimi, yayınlanmış **insan** verisine oturtulur:
Vockley ve ark. 2023 (Nat Metab, Synpheny-1 faz 2) SYNB1934 ile açlık
plazma Phe'sinde %40 düşüş bildirir. `calibrate()` bu %40'ı veren Vmax'ı
sayısal olarak bulur. Yani aktivite = 1.0 "SYNB1934 kadar iyi" demektir.

BİZİM VEKTÖR için doğru aktivite değeri BİLİNMİYOR. Suş yapılmadı,
ölçülmedi. Bu yüzden `scenarios()` bir aralık tarar; sunumda tek bir
sayı değil bu aralık gösterilmelidir.

Parametre etiketleri:
  [L] literatürden     [V] varsayım (bizim seçimimiz, ölçülmedi)
"""

from __future__ import annotations

from dataclasses import dataclass, field

# ---------------------------------------------------------------------------
# Parametreler
# ---------------------------------------------------------------------------

#: [L] Synpheny-1 faz 2 çalışmasına alınma ölçütü: tarama Phe'si > 600 µM.
#: Tedavisiz klasik PKU'da açlık Phe'si tipik olarak daha yüksektir.
#: Buradaki 1000 µM [V] bir başlangıç noktasıdır, hastaya göre değişir.
BASELINE_PHE = 1000.0

#: [L] Vockley 2023, SYNB1934: açlık plazma Phe'sinde %40 düşüş
#: (%95 GA: -52, -24). Kalibrasyon hedefi budur.
CLINICAL_DROP = 0.40
CLINICAL_DROP_CI = (0.24, 0.52)

#: [L] Aynı çalışma, SYNB1618 (ilk suş): D5-Phe AUC'de %34 düşüş;
#: SYNB1934'te %43. Adolfsen 2021: SYNB1934'ün in vivo PAL aktivitesi
#: SYNB1618'in yaklaşık iki katı.
SYNB1618_AUC_DROP = 0.34
SYNB1934_AUC_DROP = 0.43


@dataclass
class Params:
    """Modelin bütün sayıları. Hepsi tek yerde, hepsi etiketli."""
    # [V] Phe havuzunun protein sentezine geri dönüş hız sabiti.
    # Yarı ömrü ~14 saat varsayıldı; plazma Phe'si diyet değişimine
    # saatler-günler ölçeğinde cevap verir.
    k_sentez: float = 0.05
    # [V] Toplam giriş, tedavisiz denge BASELINE_PHE olacak şekilde
    # k_sentez * BASELINE_PHE'den türetilir; %60 diyet, %40 doku yıkımı.
    diyet_payi: float = 0.60
    # [V] Öğün saatleri ve emilim süresi (saat).
    ogun_saatleri: tuple[float, ...] = (8.0, 13.0, 19.0)
    emilim_suresi: float = 2.0
    # [V] Bakterinin etkin olduğu pencere: öğünden sonra kaçıncı saatte
    # başlar ve kaç saat sürer (bağırsak geçişi).
    bakteri_gecikme: float = 0.5
    bakteri_penceresi: float = 4.0
    # [V] PAL'ın Phe'ye görünür Km'i (µM). Doyum davranışını belirler.
    km: float = 1000.0
    # Kalibrasyonla bulunur (µM/sa).
    vmax: float = 0.0
    # Simülasyon adımı (saat) ve süresi (gün).
    dt: float = 0.05
    gun: int = 30

    @property
    def toplam_giris(self) -> float:
        return self.k_sentez * BASELINE_PHE

    @property
    def k_salim(self) -> float:
        return self.toplam_giris * (1.0 - self.diyet_payi)

    @property
    def ogun_basina(self) -> float:
        """Bir öğünden emilen toplam Phe (µM)."""
        return self.toplam_giris * self.diyet_payi * 24.0 / len(self.ogun_saatleri)


def _ogun_hizi(t: float, p: Params) -> float:
    """Saat t'de öğünden emilen Phe hızı (µM/sa)."""
    saat = t % 24.0
    hiz = 0.0
    for b in p.ogun_saatleri:
        if b <= saat < b + p.emilim_suresi:
            hiz += p.ogun_basina / p.emilim_suresi
    return hiz


def _bakteri_acik(t: float, p: Params) -> bool:
    """Bakteri bu anda bağırsakta Phe tüketiyor mu?"""
    saat = t % 24.0
    for b in p.ogun_saatleri:
        bas = b + p.bakteri_gecikme
        if bas <= saat < bas + p.bakteri_penceresi:
            return True
    return False


# ---------------------------------------------------------------------------
# Çözücü
# ---------------------------------------------------------------------------


@dataclass
class Run:
    t: list[float] = field(default_factory=list)      # saat
    c: list[float] = field(default_factory=list)      # µM
    aktivite: float = 0.0

    @property
    def aclik_phe(self) -> float:
        """Son tam günün sabah 7:00 değeri (ilk öğünden hemen önce).

        Klinik çalışmalarda bildirilen "açlık plazma Phe" buna karşılık
        gelir: gece boyu aç kalmış hastadan sabah alınan kan.
        """
        hedef = [i for i, t in enumerate(self.t) if abs(t % 24.0 - 7.0) < 1e-9]
        if not hedef:
            hedef = [min(range(len(self.t)),
                         key=lambda i: abs(self.t[i] % 24.0 - 7.0))]
        return self.c[hedef[-1]]

    @property
    def gunluk_aclik(self) -> list[tuple[float, float]]:
        """Her günün 07:00 değeri: [(gün, µM), ...]."""
        out = []
        for i, t in enumerate(self.t):
            if abs(t % 24.0 - 7.0) < 1e-9:
                out.append((t / 24.0, self.c[i]))
        return out

    @property
    def son_gun_ortalama(self) -> float:
        kesim = self.t[-1] - 24.0
        v = [c for t, c in zip(self.t, self.c) if t >= kesim]
        return sum(v) / len(v)

    @property
    def tepe(self) -> float:
        kesim = self.t[-1] - 24.0
        return max(c for t, c in zip(self.t, self.c) if t >= kesim)


def simulate(p: Params, aktivite: float,
             baslangic_gunu: float = 0.0) -> Run:
    """Verilen bakteri aktivitesiyle plazma Phe seyrini hesapla.

    aktivite = 1.0  ->  SYNB1934 kadar etkin (kalibrasyon noktası)
    aktivite = 0.0  ->  tedavi yok
    baslangic_gunu  ->  tedavinin başladığı gün (öncesinde tedavi yok)
    """
    r = Run(aktivite=aktivite)
    c = BASELINE_PHE
    n = int(p.gun * 24.0 / p.dt)
    t_bas = baslangic_gunu * 24.0
    for i in range(n + 1):
        t = i * p.dt
        r.t.append(t)
        r.c.append(c)
        giris = _ogun_hizi(t, p) + p.k_salim
        cikis = p.k_sentez * c
        if aktivite > 0 and t >= t_bas and _bakteri_acik(t, p):
            cikis += aktivite * p.vmax * c / (p.km + c)
        c = max(0.0, c + (giris - cikis) * p.dt)
    return r


def calibrate(p: Params | None = None,
              hedef_dusus: float = CLINICAL_DROP) -> Params:
    """Aktivite 1.0'da açlık Phe'sini `hedef_dusus` kadar düşüren Vmax'ı bul.

    İkiye bölme ile aranır; Vmax arttıkça düşüş monotondur.
    """
    p = p or Params()
    # Karşılaştırma noktası, modelin kendi tedavisiz açlık değeridir;
    # klinik çalışma da düşüşü hastanın kendi başlangıcına göre verir.
    taban = simulate(p, 0.0).aclik_phe
    hedef = taban * (1.0 - hedef_dusus)

    def aclik(vmax: float) -> float:
        q = Params(**{**p.__dict__, "vmax": vmax})
        return simulate(q, 1.0).aclik_phe

    alt, ust = 0.0, 10.0
    while aclik(ust) > hedef and ust < 1e5:
        ust *= 2.0
    for _ in range(60):
        orta = (alt + ust) / 2.0
        if aclik(orta) > hedef:
            alt = orta
        else:
            ust = orta
    return Params(**{**p.__dict__, "vmax": (alt + ust) / 2.0})


def scenarios(aktiviteler: tuple[float, ...] = (0.0, 0.25, 0.5, 1.0, 2.0),
              p: Params | None = None) -> tuple[Params, list[Run]]:
    """Kalibre edilmiş modeli bir aktivite aralığında çalıştır."""
    p = calibrate(p)
    return p, [simulate(p, a) for a in aktiviteler]


def dusus_yuzdesi(runs: list[Run]) -> dict[float, float]:
    """Her aktivite için tedavisiz duruma göre açlık Phe düşüşü (%)."""
    taban = next((r.aclik_phe for r in runs if r.aktivite == 0.0), None)
    if taban is None:
        taban = max(r.aclik_phe for r in runs)
    return {r.aktivite: (taban - r.aclik_phe) / taban * 100.0 for r in runs}
